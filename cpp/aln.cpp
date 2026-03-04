#include "aln.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <numeric>
#include <tuple>
#include <math.h>
#include "chain.hpp"
#include "clip.hpp"
#include "sv_tags.hpp"
#include "pair_scoring.hpp"
#include "revcomp.hpp"
#include "timer.hpp"
#include "nam.hpp"
#include "paf.hpp"
#include "aligner.hpp"
#include "logger.hpp"

using namespace klibpp;

static Logger& logger = Logger::get();


namespace {

/*
 * Given a primary alignment and a list of candidate alignments, select
 * supplementary alignments that cover different query regions.
 * Candidates that overlap the primary (or previously accepted supplementary)
 * query region by more than max_supp_overlap bp are rejected.
 * Candidates whose entire aligned query region is already covered are also
 * rejected.
 */
std::pair<std::vector<Alignment>, std::vector<int>> select_supplementary(
    const Alignment& primary,
    const std::vector<Alignment>& candidates,
    unsigned max_supplementary,
    int max_supp_overlap
) {
    std::vector<Alignment> supplementaries;
    std::vector<int> tied_counts;
    if (max_supplementary == 0 || candidates.empty()) return {supplementaries, tied_counts};

    // Sort candidates by score (descending)
    std::vector<size_t> indices(candidates.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
        [&](size_t a, size_t b) { return candidates[a].score > candidates[b].score; });

    // Track covered query regions (start with primary)
    std::vector<std::pair<int, int>> covered;
    covered.push_back({primary.query_start, primary.query_end});

    for (size_t idx : indices) {
        if (supplementaries.size() >= max_supplementary) break;
        const auto& cand = candidates[idx];
        int cand_len = cand.query_end - cand.query_start;
        if (cand_len <= 0) continue;

        bool too_much_overlap = false;
        for (const auto& [cov_start, cov_end] : covered) {
            int overlap = std::max(0, std::min(cand.query_end, cov_end) - std::max(cand.query_start, cov_start));
            if (overlap > max_supp_overlap || overlap >= cand_len) {
                too_much_overlap = true;
                break;
            }
        }
        if (!too_much_overlap) {
            // Count how many other candidates with the same score would also
            // pass the overlap check (i.e. are true alternative placements)
            int tied = 0;
            for (size_t j : indices) {
                if (j == idx || candidates[j].score != cand.score) continue;
                const auto& other = candidates[j];
                int other_len = other.query_end - other.query_start;
                if (other_len <= 0) continue;
                bool other_overlaps = false;
                for (const auto& [cov_start, cov_end] : covered) {
                    int ov = std::max(0, std::min(other.query_end, cov_end) - std::max(other.query_start, cov_start));
                    if (ov > max_supp_overlap || ov >= other_len) {
                        other_overlaps = true;
                        break;
                    }
                }
                if (!other_overlaps) tied++;
            }
            supplementaries.push_back(cand);
            tied_counts.push_back(tied);
            covered.push_back({cand.query_start, cand.query_end});
        }
    }
    return {supplementaries, tied_counts};
}

inline Alignment extend_seed(
    const Aligner& aligner,
    const Nam &nam,
    const References& references,
    const Read& read,
    bool consistent_nam
);

bool has_significant_soft_clip(const Alignment& primary, int read_len, int min_clip);

/*
 * Determine whether the NAM represents a match to the forward or
 * reverse-complemented sequence by checking in which orientation the
 * first and last strobe in the NAM match
 *
 * - If first and last strobe match in forward orientation, return true.
 * - If first and last strobe match in reverse orientation, update the NAM
 *   in place and return true.
 * - If first and last strobe do not match consistently, return false.
 */
bool reverse_nam_if_needed(Nam& nam, const Read& read, const References& references, int k) {
    auto read_len = read.size();
    const auto& ref_seq = references.sequences[nam.ref_id];
    const char* ref_start_ptr = ref_seq.data() + nam.ref_start;
    const char* ref_end_ptr = ref_seq.data() + nam.ref_end - k;

    const std::string& seq = nam.is_revcomp ? read.rc : read.seq;
    const std::string& seq_rc = nam.is_revcomp ? read.seq : read.rc;

    if (std::memcmp(ref_start_ptr, seq.data() + nam.query_start, k) == 0
        && std::memcmp(ref_end_ptr, seq.data() + nam.query_end - k, k) == 0) {
        return true;
    }

    // False forward or false reverse (possible due to symmetrical hash values)
    //    we need two extra checks for this - hopefully this will remove all the false hits we see (true hash collisions should be very few)
    int q_start_tmp = read_len - nam.query_end;
    int q_end_tmp = read_len - nam.query_start;
    // false reverse hit, change coordinates in nam to forward
    if (std::memcmp(ref_start_ptr, seq_rc.data() + q_start_tmp, k) == 0
        && std::memcmp(ref_end_ptr, seq_rc.data() + q_end_tmp - k, k) == 0) {
        nam.is_revcomp = !nam.is_revcomp;
        nam.query_start = q_start_tmp;
        nam.query_end = q_end_tmp;
        return true;
    }
    return false;
}

inline void align_single(
    const Aligner& aligner,
    Sam& sam,
    std::vector<Nam>& nams,
    const KSeq& record,
    int k,
    const References& references,
    Details& details,
    float dropoff_threshold,
    int max_tries,
    unsigned max_secondary,
    unsigned max_supplementary,
    int max_supp_overlap,
    std::minstd_rand& random_engine,
    bool sv_tags = false,
    SvEvidenceCollector* collector = nullptr,
    bool is_near_hotspot = false,
    const PhasingMap* phasing_map = nullptr,
    int min_clip = 15
) {
    if (nams.empty()) {
        sam.add_unmapped(record);
        return;
    }

    Read read(record.seq);
    std::vector<Alignment> alignments;
    int tries = 0;
    Nam n_max = nams[0];

    int best_edit_distance = std::numeric_limits<int>::max();
    int best_score = 0;
    int second_best_score = 0;
    int alignments_with_best_score = 0;
    size_t best_index = 0;

    Alignment best_alignment;
    best_alignment.is_unaligned = true;

    const bool collect_alignments = max_secondary > 0 || max_supplementary > 0;

    for (auto &nam : nams) {
        float score_dropoff = (float) nam.score / n_max.score;
        if (tries >= max_tries || (tries > 1 && best_edit_distance == 0) || score_dropoff < dropoff_threshold) {
            break;
        }
        bool consistent_nam = reverse_nam_if_needed(nam, read, references, k);
        details.inconsistent_nams += !consistent_nam;
        auto alignment = extend_seed(aligner, nam, references, read, consistent_nam);
        details.tried_alignment++;
        if (alignment.is_unaligned) {
            tries++;
            continue;
        }
        details.gapped += alignment.gapped;

        if (collect_alignments) {
            alignments.emplace_back(alignment);
        }

        if (alignment.score >= best_score) {
            second_best_score = best_score;
            bool update_best = false;
            if (alignment.score > best_score) {
                alignments_with_best_score = 1;
                update_best = true;
            } else {
                assert(alignment.score == best_score);
                alignments_with_best_score++;
                std::uniform_int_distribution<> distrib(1, alignments_with_best_score);
                if (distrib(random_engine) == 1) {
                    update_best = true;
                }
            }
            if (update_best) {
                best_score = alignment.score;
                best_alignment = std::move(alignment);
                best_index = collect_alignments ? alignments.size() - 1 : tries;
                best_edit_distance = best_alignment.global_ed;
            }
        } else if (alignment.score > second_best_score) {
            second_best_score = alignment.score;
        }
        tries++;
    }
    if (best_alignment.is_unaligned) {
        sam.add_unmapped(record);
        return;
    }
    details.best_alignments = alignments_with_best_score;
    uint8_t mapq = (60.0 * (best_score - second_best_score) + best_score - 1) / best_score;

    // Detect supplementary alignments before outputting primary (SA tag on
    // primary must reference all supplementaries).
    // select_supplementary() naturally skips the primary itself because its
    // query region is already in the covered set, so overlap >= cand_len.
    std::vector<Alignment> supplementaries;
    std::vector<int> supp_tied_counts;

    if (max_supplementary > 0 && alignments.size() > 1
        && has_significant_soft_clip(best_alignment, record.seq.size(), min_clip)) {
        std::tie(supplementaries, supp_tied_counts) = select_supplementary(best_alignment, alignments, max_supplementary, max_supp_overlap);
    }

    // Build SA tags if supplementaries were found
    std::string primary_sa_tag;
    std::vector<std::string> supp_sa_tags;

    if (!supplementaries.empty()) {
        std::string primary_entry = sam.format_sa_entry(best_alignment, mapq);

        // Build SA entries for each supplementary
        std::vector<std::string> supp_entries;
        for (const auto& s : supplementaries) {
            supp_entries.push_back(sam.format_sa_entry(s, mapq));
        }

        // Primary's SA tag = all supplementary entries
        primary_sa_tag = "SA:Z:";
        for (const auto& e : supp_entries) {
            primary_sa_tag.append(e);
        }

        // Each supplementary's SA tag = primary + other supplementary entries
        for (size_t i = 0; i < supplementaries.size(); ++i) {
            std::string tag = "SA:Z:";
            tag.append(primary_entry);
            for (size_t j = 0; j < supplementaries.size(); ++j) {
                if (j != i) {
                    tag.append(supp_entries[j]);
                }
            }
            supp_sa_tags.push_back(tag);
        }
    }

    // Set XS (best non-primary alignment score) from extension scores
    details.xs = second_best_score;

    // Build extra_tags for primary: SA + XD/XR + XA
    std::string extra_tags = primary_sa_tag;

    // Collect SV evidence for two-pass (pass 1) — cheap position extraction
    if (collector && !supplementaries.empty()) {
        for (const auto& s : supplementaries) {
            collector->add(best_alignment.ref_id, best_alignment.ref_start + best_alignment.length);
            collector->add(s.ref_id, s.ref_start);
        }
    }

    // SV tags: XD/XR/XE/XF/XU + YS + breakpoint tags (--sv-tags only, skip in pass 1)
    std::vector<BreakpointInfo> bp_infos;
    if (sv_tags && !supplementaries.empty() && !collector) {
        // XD/XR: score delta and ratio between primary and best supplementary
        float score_ratio = 0.0f;
        {
            auto dt = compute_supp_delta_tags(best_alignment, supplementaries);
            if (!dt.empty()) {
                if (!extra_tags.empty()) extra_tags += '\t';
                extra_tags += dt;
            }
            if (best_alignment.score > 0 && !supplementaries.empty()) {
                score_ratio = static_cast<float>(supplementaries[0].score) / best_alignment.score;
            }
        }
        // XE: split chain entropy over all alignments
        float xe = compute_split_entropy(alignments);
        {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "XE:f:%.3f", xe);
            if (!extra_tags.empty()) extra_tags += '\t';
            extra_tags += buf;
        }
        // XF: max artifact probability across supplementaries
        float xf = compute_max_artifact_prob(best_alignment, supplementaries, record.seq.length());
        {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "XF:f:%.2f", xf);
            if (!extra_tags.empty()) extra_tags += '\t';
            extra_tags += buf;
        }
        // Compute breakpoint info for each supplementary
        for (const auto& s : supplementaries) {
            bp_infos.push_back(compute_breakpoint_info(best_alignment, s, references, record.seq));
        }
        // YS: comma-separated SV types on primary
        {
            if (!extra_tags.empty()) extra_tags += '\t';
            extra_tags += "YS:Z:";
            for (size_t i = 0; i < bp_infos.size(); ++i) {
                if (i > 0) extra_tags += ',';
                extra_tags += bp_infos[i].sv_type;
            }
        }
        // Breakpoint tags from best supplementary on primary
        if (!bp_infos.empty()) {
            extra_tags += '\t';
            extra_tags += bp_infos[0].format_all(references);
        }
        // XU: REF/SV/RPT posteriors
        {
            int max_tied = 0;
            for (int tc : supp_tied_counts) {
                if (tc > max_tied) max_tied = tc;
            }
            float total_hits = static_cast<float>(details.hits.total_hits());
            float filtered_hits = static_cast<float>(details.hits.total_filtered());
            float mappability = total_hits > 0 ? 1.0f - (filtered_hits / total_hits) : 1.0f;
            auto xu = compute_xu_posteriors(score_ratio, xf, xe, mappability, max_tied);
            if (!extra_tags.empty()) extra_tags += '\t';
            extra_tags += format_xu_tag(xu);
        }
    }

    // Tag hotspot-realigned reads
    if (is_near_hotspot) {
        if (!extra_tags.empty()) extra_tags += '\t';
        extra_tags += "XH:i:1";
    }

    // XA: alternative alignment locations (all non-primary, non-supplementary, deduplicated)
    if (collect_alignments && alignments.size() > 1) {
        std::string xa_entries;
        // Track seen (ref_id, ref_start, strand) to avoid duplicates
        std::vector<std::tuple<int, int, bool>> seen_xa;
        for (size_t i = 0; i < alignments.size(); ++i) {
            if (i == best_index) continue;
            const auto& aln = alignments[i];
            bool is_supp = false;
            for (const auto& s : supplementaries) {
                if (aln.ref_id == s.ref_id && aln.ref_start == s.ref_start
                    && aln.query_start == s.query_start && aln.query_end == s.query_end) {
                    is_supp = true;
                    break;
                }
            }
            if (is_supp) continue;
            auto key = std::make_tuple(aln.ref_id, aln.ref_start, aln.is_revcomp);
            bool dup = false;
            for (const auto& k : seen_xa) {
                if (k == key) { dup = true; break; }
            }
            if (dup) continue;
            seen_xa.push_back(key);
            xa_entries += sam.format_xa_entry(aln);
        }
        if (!xa_entries.empty()) {
            if (!extra_tags.empty()) extra_tags += '\t';
            extra_tags += "XA:Z:";
            extra_tags += xa_entries;
        }
    }

    // XK: soft-clip realignment for short clips not covered by supplementaries
    if (sv_tags && !collector) {
        bool left_covered = false, right_covered = false;
        for (const auto& s : supplementaries) {
            if (s.query_start < best_alignment.query_start) left_covered = true;
            if (s.query_end > best_alignment.query_end) right_covered = true;
        }
        // Coverage flags are in forward-strand coords; for revcomp reads,
        // forward-strand left = mapped-strand trailing (3') and vice versa.
        if (best_alignment.is_revcomp) std::swap(left_covered, right_covered);
        const auto& mapped_seq = best_alignment.is_revcomp ? read.rc : read.seq;
        auto xk_clips = realign_soft_clips(
            best_alignment, mapped_seq, references, aligner,
            left_covered, right_covered);
        if (!xk_clips.empty()) {
            if (!extra_tags.empty()) extra_tags += '\t';
            extra_tags += "XK:Z:" + format_xk_tag(xk_clips);
        }
    }

    // HP/PS: haplotype phasing
    if (phasing_map && !best_alignment.is_unaligned) {
        auto phase = phasing_map->assign_haplotype(best_alignment, record.seq, best_alignment.cigar);
        if (phase) {
            if (!extra_tags.empty()) extra_tags += '\t';
            extra_tags += format_phase_tags(*phase);
        }
    }

    // Output primary alignment
    sam.add(best_alignment, record, read.rc, mapq, PRIMARY, details, extra_tags);

    // Output supplementary alignments
    for (size_t i = 0; i < supplementaries.size(); ++i) {
        std::string supp_extra = supp_sa_tags[i];
        if (sv_tags && i < bp_infos.size()) {
            if (!supp_extra.empty()) supp_extra += '\t';
            supp_extra += "YS:Z:" + bp_infos[i].sv_type;
            supp_extra += '\t';
            supp_extra += bp_infos[i].format_all(references);
        }
        if (sv_tags && i < supp_tied_counts.size() && supp_tied_counts[i] > 0) {
            if (!supp_extra.empty()) supp_extra += '\t';
            supp_extra += "Xt:i:" + std::to_string(supp_tied_counts[i]);
        }
        sam.add(supplementaries[i], record, read.rc, mapq, SUPPLEMENTARY_ALN, details, supp_extra);
    }

    if (max_secondary == 0 || alignments.empty()) {
        return;
    }

    // Secondary alignments

    // Remove the primary alignment
    if (alignments.size() > 1) {
        std::swap(alignments[best_index], alignments[alignments.size() - 1]);
    }
    alignments.resize(alignments.size() - 1);

    // Sort remaining alignments by score, highest first
    std::sort(alignments.begin(), alignments.end(),
        [](const Alignment& a, const Alignment& b) { return a.score > b.score; });

    // Output secondary alignments (skip those already output as supplementary)
    double secondary_dropoff = 2 * aligner.parameters.mismatch + aligner.parameters.gap_open;
    size_t n = 0;
    for (const auto& aln : alignments) {
        if (n >= max_secondary || best_score - aln.score > secondary_dropoff) {
            break;
        }
        // Skip alignments already emitted as supplementary
        bool is_supp = false;
        for (const auto& s : supplementaries) {
            if (aln.ref_id == s.ref_id && aln.ref_start == s.ref_start
                && aln.query_start == s.query_start && aln.query_end == s.query_end) {
                is_supp = true;
                break;
            }
        }
        if (is_supp) continue;
        sam.add(aln, record, read.rc, mapq, SECONDARY_ALN, details);
        n++;
    }
}

/*
 Extend a NAM so that it covers the entire read and return the resulting
 alignment.
*/
inline Alignment extend_seed(
    const Aligner& aligner,
    const Nam &nam,
    const References& references,
    const Read& read,
    bool consistent_nam
) {
    const std::string& query = nam.is_revcomp ? read.rc : read.seq;
    const std::string& ref = references.sequences[nam.ref_id];

    const auto projected_ref_start = nam.projected_ref_start();
    const auto projected_ref_end = std::min(nam.ref_end + query.size() - nam.query_end, ref.size());

    AlignmentInfo info;
    int result_ref_start;
    bool gapped = true;
    if (projected_ref_end - projected_ref_start == query.size() && consistent_nam) {
        std::string ref_segm_ham = ref.substr(projected_ref_start, query.size());
        auto hamming_dist = hamming_distance(query, ref_segm_ham);

        if (hamming_dist >= 0 && (((float) hamming_dist / query.size()) < 0.05) ) { //Hamming distance worked fine, no need to ksw align
            info = hamming_align(query, ref_segm_ham, aligner.parameters.match, aligner.parameters.mismatch, aligner.parameters.end_bonus);
            result_ref_start = projected_ref_start + info.ref_start;
            gapped = false;
        }
    }
    if (gapped) {
        const int diff = std::abs(nam.ref_span() - nam.query_span());
        const int ext_left = std::min(50, projected_ref_start);
        const int ref_start = projected_ref_start - ext_left;
        const int ext_right = std::min(std::size_t(50), ref.size() - nam.ref_end);
        const auto ref_segm_size = read.size() + diff + ext_left + ext_right;
        const auto ref_segm = ref.substr(ref_start, ref_segm_size);
        auto opt_info = aligner.align(query, ref_segm);
        if (opt_info) {
            info = opt_info.value();
            result_ref_start = ref_start + info.ref_start;
        } else {
            // TODO This function should instead return an std::optional<Alignment>
            Alignment alignment;
            alignment.is_unaligned = true;
            alignment.edit_distance = 100000;
            alignment.ref_start = 0;
            alignment.score = -100000;

            return alignment;
        }
    }
    int softclipped = info.query_start + (query.size() - info.query_end);
    Alignment alignment;
    alignment.cigar = std::move(info.cigar);
    alignment.edit_distance = info.edit_distance;
    alignment.global_ed = info.edit_distance + softclipped;
    alignment.score = info.sw_score;
    alignment.ref_start = result_ref_start;
    alignment.length = info.ref_span();
    alignment.is_revcomp = nam.is_revcomp;
    alignment.is_unaligned = false;
    alignment.ref_id = nam.ref_id;
    alignment.gapped = gapped;

    // Convert query coordinates to forward-strand
    if (nam.is_revcomp) {
        int read_len = static_cast<int>(query.size());
        alignment.query_start = read_len - info.query_end;
        alignment.query_end = read_len - info.query_start;
    } else {
        alignment.query_start = info.query_start;
        alignment.query_end = info.query_end;
    }

    return alignment;
}

struct ReadSupplementaryInfo {
    std::vector<Alignment> supps;
    std::vector<int> tied_counts;  // Per-supp: how many candidates shared the same score
    std::string primary_sa_tag;  // SA tag for the primary record
    std::vector<std::string> supp_sa_tags;  // SA tag for each supplementary record
};

/*
 * Build SA tags linking primary and supplementary records.
 * Populates primary_sa_tag and supp_sa_tags in info.
 */
void build_sa_tags(
    ReadSupplementaryInfo& info,
    Sam& sam,
    const Alignment& primary,
    uint8_t mapq
) {
    if (info.supps.empty()) return;

    std::string primary_entry = sam.format_sa_entry(primary, mapq);
    std::vector<std::string> supp_entries;
    for (auto& s : info.supps) {
        supp_entries.push_back(sam.format_sa_entry(s, mapq));
    }

    info.primary_sa_tag = "SA:Z:";
    for (auto& e : supp_entries) {
        info.primary_sa_tag.append(e);
    }

    for (size_t i = 0; i < info.supps.size(); ++i) {
        std::string tag = "SA:Z:";
        tag.append(primary_entry);
        for (size_t j = 0; j < info.supps.size(); ++j) {
            if (j != i) tag.append(supp_entries[j]);
        }
        info.supp_sa_tags.push_back(tag);
    }
}

/*
 * Check if primary alignment has enough soft clipping to warrant
 * searching for supplementary alignments.
 */
bool has_significant_soft_clip(const Alignment& primary, int read_len, int min_clip) {
    int left_clip = primary.query_start;
    int right_clip = read_len - primary.query_end;
    return left_clip >= min_clip || right_clip >= min_clip;
}

/*
 * Detect supplementary alignments for one read by aligning its NAMs,
 * then build SA tags linking primary and supplementary records.
 * Used on the fast path where only the primary has been aligned.
 */
ReadSupplementaryInfo build_supplementary_info(
    Sam& sam,
    const Aligner& aligner,
    std::vector<Nam>& nams,
    const Read& read,
    const Alignment& primary,
    uint8_t mapq,
    int k,
    const References& references,
    float dropoff_threshold,
    unsigned max_tries,
    unsigned max_supplementary,
    int max_supp_overlap,
    int min_clip
) {
    ReadSupplementaryInfo info;
    if (max_supplementary == 0 || nams.empty() || primary.is_unaligned) {
        return info;
    }

    // Early exit: if primary covers most of the read (both soft clips < min_clip),
    // there's no room for a meaningful supplementary alignment
    if (!has_significant_soft_clip(primary, read.size(), min_clip)) {
        return info;
    }

    // Build candidate alignments from NAMs
    std::vector<Alignment> candidates;
    auto n_max_score = nams[0].score;
    unsigned tries = 0;
    for (auto& nam : nams) {
        if (tries >= max_tries) break;
        float score_dropoff = (float)nam.score / n_max_score;
        if (tries > 1 && score_dropoff < dropoff_threshold) break;

        bool consistent = reverse_nam_if_needed(nam, read, references, k);
        auto aln = extend_seed(aligner, nam, references, read, consistent);
        if (!aln.is_unaligned) {
            candidates.push_back(aln);
        }
        tries++;
    }

    std::tie(info.supps, info.tied_counts) = select_supplementary(primary, candidates, max_supplementary, max_supp_overlap);
    build_sa_tags(info, sam, primary, mapq);

    return info;
}

/*
 * Build supplementary info from pre-computed candidate alignments.
 * Used on the normal paired-end path where alignments are already available.
 */
ReadSupplementaryInfo build_supplementary_info_from_candidates(
    Sam& sam,
    const std::vector<Alignment>& candidates,
    const Alignment& primary,
    uint8_t mapq,
    unsigned max_supplementary,
    int max_supp_overlap,
    int read_len,
    int min_clip
) {
    ReadSupplementaryInfo info;
    if (max_supplementary == 0 || primary.is_unaligned || candidates.empty()) {
        return info;
    }

    // Early exit: if primary covers most of the read, skip
    if (!has_significant_soft_clip(primary, read_len, min_clip)) {
        return info;
    }

    std::tie(info.supps, info.tied_counts) = select_supplementary(primary, candidates, max_supplementary, max_supp_overlap);
    build_sa_tags(info, sam, primary, mapq);

    return info;
}


/*
 * Align a pair of reads for which only one has NAMs. For the other, rescue
 * is attempted by aligning it locally.
 */
std::vector<ScoredAlignmentPair> rescue_read(
    const Read& read2,  // read to be rescued
    const Read& read1,  // read that has NAMs
    const Aligner& aligner,
    const References& references,
    std::vector<Nam> &nams1,
    int max_tries,
    float dropoff,
    std::array<Details, 2>& details,
    int k,
    float mu,
    float sigma,
    float rescue_sigma_mult
) {
    Nam n_max1 = nams1[0];
    int tries = 0;

    std::vector<Alignment> alignments1;
    std::vector<Alignment> alignments2;
    for (auto& nam : nams1) {
        float score_dropoff1 = (float) nam.n_matches / n_max1.n_matches;
        // only consider top hits (as minimap2 does) and break if below dropoff cutoff.
        if (tries >= max_tries || score_dropoff1 < dropoff) {
            break;
        }

        const bool consistent_nam = reverse_nam_if_needed(nam, read1, references, k);
        details[0].inconsistent_nams += !consistent_nam;
        auto alignment = extend_seed(aligner, nam, references, read1, consistent_nam);
        details[0].gapped += alignment.gapped;
        alignments1.emplace_back(alignment);
        details[0].tried_alignment++;

        // Force SW alignment to rescue mate
        Alignment a2 = rescue_align(aligner, nam, references, read2, mu, sigma, k, rescue_sigma_mult);
        details[1].mate_rescue += !a2.is_unaligned;
        alignments2.emplace_back(a2);

        tries++;
    }
    std::sort(alignments1.begin(), alignments1.end(), by_score<Alignment>);
    std::sort(alignments2.begin(), alignments2.end(), by_score<Alignment>);

    // Calculate best combined score here
    auto high_scores = get_best_scoring_pairs(alignments1, alignments2, mu, sigma );

    return high_scores;
}

void output_aligned_pairs(
    const std::vector<ScoredAlignmentPair>& high_scores,
    Sam& sam,
    size_t max_secondary,
    double secondary_dropoff,
    const KSeq& record1,
    const KSeq& record2,
    const Read& read1,
    const Read& read2,
    float mu,
    float sigma,
    const std::array<Details, 2>& details,
    bool sv_tags,
    const PhasingMap* phasing_map = nullptr
) {

    if (high_scores.empty()) {
        sam.add_unmapped_pair(record1, record2);
        return;
    }

    auto [mapq1, mapq2] = joint_mapq_from_high_scores(high_scores);
    auto best_aln_pair = high_scores[0];

    // append both alignments to string here
    if (max_secondary == 0) {
        Alignment alignment1 = best_aln_pair.alignment1;
        Alignment alignment2 = best_aln_pair.alignment2;

        std::string extra1, extra2;
        if (sv_tags) {
            extra1 = build_sv_tags(alignment1, alignment2, mu, sigma, details[0], false);
            extra2 = build_sv_tags(alignment1, alignment2, mu, sigma, details[1], false);
        }
        if (phasing_map) {
            if (!alignment1.is_unaligned) {
                auto p1 = phasing_map->assign_haplotype(alignment1, record1.seq, alignment1.cigar);
                if (p1) { if (!extra1.empty()) extra1 += '\t'; extra1 += format_phase_tags(*p1); }
            }
            if (!alignment2.is_unaligned) {
                auto p2 = phasing_map->assign_haplotype(alignment2, record2.seq, alignment2.cigar);
                if (p2) { if (!extra2.empty()) extra2 += '\t'; extra2 += format_phase_tags(*p2); }
            }
        }
        sam.add_pair(alignment1, alignment2, record1, record2, read1.rc, read2.rc, mapq1, mapq2, is_proper_pair(alignment1, alignment2, mu, sigma), PRIMARY, details, extra1, extra2);
    } else {
        auto max_out = std::min(high_scores.size(), max_secondary);
        AlignmentType aln_type = PRIMARY;
        float s_max = best_aln_pair.score;
        for (size_t i = 0; i < max_out; ++i) {
            auto aln_pair = high_scores[i];
            Alignment alignment1 = aln_pair.alignment1;
            Alignment alignment2 = aln_pair.alignment2;
            float s_score = aln_pair.score;
            if (i > 0) {
                aln_type = SECONDARY_ALN;
                mapq1 = 0;
                mapq2 = 0;
            }
            if (s_max - s_score < secondary_dropoff) {
                bool is_proper = is_proper_pair(alignment1, alignment2, mu, sigma);
                std::string extra1, extra2;
                if (sv_tags && i == 0) {
                    extra1 = build_sv_tags(alignment1, alignment2, mu, sigma, details[0], false);
                    extra2 = build_sv_tags(alignment1, alignment2, mu, sigma, details[1], false);
                }
                if (phasing_map && i == 0) {
                    if (!alignment1.is_unaligned) {
                        auto p1 = phasing_map->assign_haplotype(alignment1, record1.seq, alignment1.cigar);
                        if (p1) { if (!extra1.empty()) extra1 += '\t'; extra1 += format_phase_tags(*p1); }
                    }
                    if (!alignment2.is_unaligned) {
                        auto p2 = phasing_map->assign_haplotype(alignment2, record2.seq, alignment2.cigar);
                        if (p2) { if (!extra2.empty()) extra2 += '\t'; extra2 += format_phase_tags(*p2); }
                    }
                }
                sam.add_pair(alignment1, alignment2, record1, record2, read1.rc, read2.rc, mapq1, mapq2, is_proper, aln_type, details, extra1, extra2);
            } else {
                break;
            }
        }
    }
}

// compute dropoff of the first (top) NAM
float top_dropoff(std::vector<Nam>& nams) {
    auto& n_max = nams[0];
    if (n_max.n_matches <= 2) {
        return 1.0;
    }
    if (nams.size() > 1) {
        return (float) nams[1].n_matches / n_max.n_matches;
    }
    return 0.0;
}

std::vector<ScoredAlignmentPair> align_paired(
    const Aligner& aligner,
    std::vector<Nam> &nams1,
    std::vector<Nam> &nams2,
    const Read& read1,
    const Read& read2,
    int k,
    const References& references,
    std::array<Details, 2>& details,
    float dropoff,
    const InsertSizeDistribution &isize_est,
    unsigned max_tries,
    float rescue_sigma_mult = 5.0f
) {
    const auto mu = isize_est.mu;
    const auto sigma = isize_est.sigma;

    if (nams1.empty() && nams2.empty()) {
         // None of the reads have any NAMs
        return std::vector<ScoredAlignmentPair>{};
    }

    if (!nams1.empty() && nams2.empty()) {
        // Only read 1 has NAMS: attempt to rescue read 2
        return rescue_read(
            read2,
            read1,
            aligner,
            references,
            nams1,
            max_tries,
            dropoff,
            details,
            k,
            mu,
            sigma,
            rescue_sigma_mult
        );
    }

    if (nams1.empty() && !nams2.empty()) {
        // Only read 2 has NAMS: attempt to rescue read 1
        std::array<Details, 2> swapped_details{details[1], details[0]};
        std::vector<ScoredAlignmentPair> pairs = rescue_read(
            read1,
            read2,
            aligner,
            references,
            nams2,
            max_tries,
            dropoff,
            swapped_details,
            k,
            mu,
            sigma,
            rescue_sigma_mult
        );
        details[0] += swapped_details[1];
        details[1] += swapped_details[0];
        for (auto& pair : pairs) {
            std::swap(pair.alignment1, pair.alignment2);
        }

        return pairs;
    }

    // If we get here, both reads have NAMs
    assert(!nams1.empty() && !nams2.empty());

    // Deal with the typical case that both reads map uniquely and form a proper pair
    if (top_dropoff(nams1) < dropoff && top_dropoff(nams2) < dropoff && is_proper_nam_pair(nams1[0], nams2[0], mu, sigma)) {
        Nam n_max1 = nams1[0];
        Nam n_max2 = nams2[0];

        bool consistent_nam1 = reverse_nam_if_needed(n_max1, read1, references, k);
        details[0].inconsistent_nams += !consistent_nam1;
        bool consistent_nam2 = reverse_nam_if_needed(n_max2, read2, references, k);
        details[1].inconsistent_nams += !consistent_nam2;

        auto alignment1 = extend_seed(aligner, n_max1, references, read1, consistent_nam1);
        details[0].tried_alignment++;
        details[0].gapped += alignment1.gapped;
        auto alignment2 = extend_seed(aligner, n_max2, references, read2, consistent_nam2);
        details[1].tried_alignment++;
        details[1].gapped += alignment2.gapped;

        return std::vector<ScoredAlignmentPair>{{-1, alignment1, alignment2}};
    }

    // Do a full search for highest-scoring pair
    // Get top hit counts for all locations. The joint hit count is the sum of hits of the two mates. Then align as long as score dropoff or cnt < 20

    std::vector<NamPair> nam_pairs = get_best_scoring_nam_pairs(nams1, nams2, mu, sigma);

    // Cache for already computed alignments. Maps NAM ids to alignments.
    robin_hood::unordered_map<int,Alignment> is_aligned1;
    robin_hood::unordered_map<int,Alignment> is_aligned2;

    // These keep track of the alignments that would be best if we treated
    // the paired-end read as two single-end reads.
    Alignment a1_indv_max, a2_indv_max;
    {
        auto n1_max = nams1[0];
        bool consistent_nam1 = reverse_nam_if_needed(n1_max, read1, references, k);
        details[0].inconsistent_nams += !consistent_nam1;
        a1_indv_max = extend_seed(aligner, n1_max, references, read1, consistent_nam1);
        is_aligned1[n1_max.nam_id] = a1_indv_max;
        details[0].tried_alignment++;
        details[0].gapped += a1_indv_max.gapped;

        auto n2_max = nams2[0];
        bool consistent_nam2 = reverse_nam_if_needed(n2_max, read2, references, k);
        details[1].inconsistent_nams += !consistent_nam2;
        a2_indv_max = extend_seed(aligner, n2_max, references, read2, consistent_nam2);
        is_aligned2[n2_max.nam_id] = a2_indv_max;
        details[1].tried_alignment++;
        details[1].gapped += a2_indv_max.gapped;
    }

    // Turn pairs of high-scoring NAMs into pairs of alignments
    std::vector<ScoredAlignmentPair> high_scores;
    auto max_score = nam_pairs[0].score;
    for (auto &[score_, n1, n2] : nam_pairs) {
        float score_dropoff = (float) score_ / max_score;

        if (high_scores.size() >= max_tries || score_dropoff < dropoff) {
            break;
        }

        // Get alignments for the two NAMs, either by computing the alignment,
        // retrieving it from the cache or by attempting a rescue (if the NAM
        // actually is a dummy, that is, only the partner is available)
        Alignment a1;
        // ref_start == -1 is a marker for a dummy NAM
        if (n1.ref_start >= 0) {
            if (is_aligned1.find(n1.nam_id) != is_aligned1.end() ){
                a1 = is_aligned1[n1.nam_id];
            } else {
                bool consistent_nam = reverse_nam_if_needed(n1, read1, references, k);
                details[0].inconsistent_nams += !consistent_nam;
                a1 = extend_seed(aligner, n1, references, read1, consistent_nam);
                is_aligned1[n1.nam_id] = a1;
                details[0].tried_alignment++;
                details[0].gapped += a1.gapped;
            }
        } else {
            details[1].inconsistent_nams += !reverse_nam_if_needed(n2, read2, references, k);
            a1 = rescue_align(aligner, n2, references, read1, mu, sigma, k, rescue_sigma_mult);
            details[0].mate_rescue += !a1.is_unaligned;
            details[0].tried_alignment++;
        }
        if (a1.score > a1_indv_max.score) {
            a1_indv_max = a1;
        }

        Alignment a2;
        // ref_start == -1 is a marker for a dummy NAM
        if (n2.ref_start >= 0) {
            if (is_aligned2.find(n2.nam_id) != is_aligned2.end() ){
                a2 = is_aligned2[n2.nam_id];
            } else {
                bool consistent_nam = reverse_nam_if_needed(n2, read2, references, k);
                details[1].inconsistent_nams += !consistent_nam;
                a2 = extend_seed(aligner, n2, references, read2, consistent_nam);
                is_aligned2[n2.nam_id] = a2;
                details[1].tried_alignment++;
                details[1].gapped += a2.gapped;
            }
        } else {
            details[0].inconsistent_nams += !reverse_nam_if_needed(n1, read1, references, k);
            a2 = rescue_align(aligner, n1, references, read2, mu, sigma, k, rescue_sigma_mult);
            details[1].mate_rescue += !a2.is_unaligned;
            details[1].tried_alignment++;
        }
        if (a2.score > a2_indv_max.score){
            a2_indv_max = a2;
        }

        bool r1_r2 = a2.is_revcomp && (a1.ref_start <= a2.ref_start) && ((a2.ref_start - a1.ref_start) < mu + 10*sigma); // r1 ---> <---- r2
        bool r2_r1 = a1.is_revcomp && (a2.ref_start <= a1.ref_start) && ((a1.ref_start - a2.ref_start) < mu + 10*sigma); // r2 ---> <---- r1

        double combined_score;
        if (r1_r2 || r2_r1) {
            // Treat a1/a2 as a pair
            float x = std::abs(a1.ref_start - a2.ref_start);
            combined_score = (double)a1.score + (double)a2.score + std::max(-20.0f + 0.001f, log(normal_pdf(x, mu, sigma)));
            //* (1 - s2 / s1) * min_matches * log(s1);
        } else {
            // Treat a1/a2 as two single-end reads
            // 20 corresponds to a value of log(normal_pdf(x, mu, sigma)) of more than 5 stddevs away (for most reasonable values of stddev)
            combined_score = (double)a1.score + (double)a2.score - 20;
        }

        ScoredAlignmentPair aln_pair{combined_score, a1, a2};
        high_scores.push_back(aln_pair);
    }

    // Finally, add highest scores of both mates as individually mapped
    double combined_score = (double)a1_indv_max.score + (double)a2_indv_max.score - 20; // 20 corresponds to  a value of log( normal_pdf(x, mu, sigma ) ) of more than 5 stddevs away (for most reasonable values of stddev)
    ScoredAlignmentPair aln_tuple{combined_score, a1_indv_max, a2_indv_max};
    high_scores.push_back(aln_tuple);

    return high_scores;
}

// Used for PAF and abundances output
inline void get_best_map_location(
    std::vector<Nam> &nams1,
    std::vector<Nam> &nams2,
    InsertSizeDistribution &isize_est,
    Nam &best_nam1,
    Nam &best_nam2,
    int read1_len,
    int read2_len,
    std::vector<double> &abundances,
    bool output_abundance
) {
    std::vector<NamPair> nam_pairs = get_best_scoring_nam_pairs(nams1, nams2, isize_est.mu, isize_est.sigma);
    best_nam1.ref_start = -1; //Unmapped until proven mapped
    best_nam2.ref_start = -1; //Unmapped until proven mapped

    if (nam_pairs.empty()) {
        return;
    }

    // get best joint score
    float score_joint = 0;
    Nam n1_joint_max, n2_joint_max;
    for (auto &[score, nam1, nam2] : nam_pairs) { // already sorted by descending score
        if (nam1.ref_start >= 0 && nam2.ref_start >=0) { // Valid pair
            score_joint = nam1.score + nam2.score;
            n1_joint_max = nam1;
            n2_joint_max = nam2;
            break;
        }
    }

    // get individual best scores
    float score_indiv = 0;
    if (!nams1.empty()) {
        score_indiv += nams1[0].score / 2.0; //Penalty for being mapped individually
        best_nam1 = nams1[0];
    }
    if (!nams2.empty()) {
        score_indiv += nams2[0].score / 2.0; //Penalty for being mapped individually
        best_nam2 = nams2[0];
    }
    if (score_joint > score_indiv) { // joint score is better than individual
        best_nam1 = n1_joint_max;
        best_nam2 = n2_joint_max;

        if (output_abundance){
            // we loop twice because we need to count the number of best pairs
            size_t n_best = 0;
            for (auto &[score, n1, n2] : nam_pairs){
                if ((n1.score + n2.score) == score_joint){
                    ++n_best;
                } else {
                    break;
                }
            }
            for (auto &[score, n1, n2] : nam_pairs){
                if ((n1.score + n2.score) == score_joint){
                    if (n1.ref_start >= 0) {
                        abundances[n1.ref_id] += float(read1_len) / float(n_best);
                    }
                    if (n2.ref_start >= 0) {
                        abundances[n2.ref_id] += float(read2_len) / float(n_best);
                    }
                } else {
                    break;
                }
            }
        }
    } else if (output_abundance) {
        for (auto &[nams, read_len]: {  std::make_pair(std::cref(nams1), read1_len),
                                        std::make_pair(std::cref(nams2), read2_len) }) {
            size_t best_score = 0;
            // We loop twice because we need to count the number of NAMs with best score
            for (auto &nam : nams) {
                if (nam.score == nams[0].score){
                    ++best_score;
                } else {
                    break;
                }
            }
            for (auto &nam: nams) {
                if (nam.ref_start < 0) {
                    continue;
                }
                if (nam.score != nams[0].score){
                    break;
                }
                abundances[nam.ref_id] += float(read_len) / float(best_score);
            }
        }
    }

    if (isize_est.sample_size < 400 && score_joint > score_indiv) {
        isize_est.update(std::abs(n1_joint_max.ref_start - n2_joint_max.ref_start));
    }
}

} // end of anonymous namespace

template <typename T>
bool by_score(const T& a, const T& b)
{
    return a.score > b.score;
}

/* Shuffle the top-scoring NAMs. Input must be sorted by score.
 *
 * This helps to ensure we pick a random location in case there are multiple
 * equally good ones.
 */
void shuffle_top_nams(std::vector<Nam>& nams, std::minstd_rand& random_engine) {
    if (nams.empty()) {
        return;
    }
    auto best_score = nams[0].score;
    auto it = std::find_if(nams.begin(), nams.end(), [&](const Nam& nam) { return nam.score != best_score; });
    if (it > nams.begin() + 1) {
        std::shuffle(nams.begin(), it, random_engine);
    }
}

/*
 * Determine (roughly) whether the read sequence has some l-mer (with l = k*2/3)
 * in common with the reference sequence
 */
bool has_shared_substring(const std::string& read_seq, const std::string& ref_seq, int k) {
    size_t sub_size = 2 * k / 3;
    size_t step_size = k / 3;
    std::string_view ref_view(ref_seq);
    std::string_view read_view(read_seq);
    for (size_t i = 0; i + sub_size <= read_view.size(); i += step_size) {
        if (ref_view.find(read_view.substr(i, sub_size)) != std::string_view::npos) {
            return true;
        }
    }
    return false;
}

std::vector<Nam> get_nams_or_chains(
    const KSeq& record,
    const StrobemerIndex& index,
    const Chainer& chainer,
    AlignmentStatistics& statistics,
    Details& details,
    const MappingParameters &map_param,
    const IndexParameters& index_parameters,
    std::minstd_rand& random_engine
) {

    // Compute randstrobes
    Timer strobe_timer;
    auto query_randstrobes = randstrobes_query(record.seq, index_parameters);
    statistics.n_randstrobes += query_randstrobes[0].size() + query_randstrobes[1].size();
    details.n_seeds = query_randstrobes[0].size() + query_randstrobes[1].size();
    statistics.tot_construct_strobemers += strobe_timer.duration();

    std::vector<Nam> nams;
    if (map_param.use_nams) {
        nams = get_nams(query_randstrobes, index, statistics, details, map_param);
    } else {
        nams = chainer.get_chains(query_randstrobes, index, statistics, details, map_param);
    }

    // Sort by score
    Timer nam_sort_timer;
    std::sort(nams.begin(), nams.end(), by_score<Nam>);
    shuffle_top_nams(nams, random_engine);
    statistics.tot_sort_nams += nam_sort_timer.duration();

    if (logger.level() <= LOG_TRACE) {
        logger.trace() << "Found " << nams.size() << (map_param.use_nams ? " NAMs\n" : " chains\n");
        uint printed = 0;
        for (const auto& nam : nams) {
            if (nam.n_matches > 1 || printed < 10) {
                logger.trace() << "- " << nam << '\n';
                printed++;
            }
        }
        if (printed != nams.size()) {
            logger.trace() << "+" << nams.size() - printed << " single-anchor chains)\n";
        }
    }

    return nams;
}

void align_or_map_paired(
    const KSeq &record1,
    const KSeq &record2,
    Sam& sam,
    std::string& outstring,
    AlignmentStatistics &statistics,
    InsertSizeDistribution &isize_est,
    const Aligner &aligner,
    const Chainer& chainer,
    const MappingParameters &map_param,
    const IndexParameters& index_parameters,
    const References& references,
    const StrobemerIndex& index,
    std::minstd_rand& random_engine,
    std::vector<double> &abundances,
    SvEvidenceCollector* collector,
    const HotspotMap* hotspot_map,
    const MappingParameters* relaxed_params,
    const PhasingMap* phasing_map
) {
    std::array<Details, 2> details;
    std::array<std::vector<Nam>, 2> nams_pair;

    for (size_t is_r1 : {0, 1}) {
        const auto& record = is_r1 == 0 ? record1 : record2;
        logger.trace() << "\nQuery: " << record.name << " (R" << is_r1 + 1 << ")\n";
        nams_pair[is_r1] = get_nams_or_chains(
            record, index, chainer, statistics, details[is_r1], map_param, index_parameters, random_engine
        );
        if (!nams_pair[is_r1].empty()) {
            details[is_r1].s1 = std::lroundf(nams_pair[is_r1][0].score);
            if (nams_pair[is_r1].size() > 1) details[is_r1].s2 = std::lroundf(nams_pair[is_r1][1].score);
            details[is_r1].cm = nams_pair[is_r1][0].n_matches;
            details[is_r1].xc = std::lroundf(nams_pair[is_r1][0].score);
        }
    }

    // Two-pass: check if either read's top NAM is near an SV hotspot
    bool is_near_hotspot = false;
    if (hotspot_map) {
        for (const auto& nams : nams_pair) {
            if (!nams.empty() && hotspot_map->near_hotspot(nams[0].ref_id, nams[0].ref_start)) {
                is_near_hotspot = true;
                break;
            }
        }
    }
    const auto& effective_params = (is_near_hotspot && relaxed_params) ? *relaxed_params : map_param;

    Timer extend_timer;
    if (map_param.output_format != OutputFormat::SAM) { // PAF or abundance
        Nam nam_read1;
        Nam nam_read2;
        get_best_map_location(
                nams_pair[0], nams_pair[1],
                isize_est,
                nam_read1, nam_read2,
                record1.seq.length(), record2.seq.length(),
                abundances,
                map_param.output_format == OutputFormat::Abundance);
        if (map_param.output_format == OutputFormat::PAF) {
            uint8_t mapq1 = proper_pair_mapq(nams_pair[0]);
            uint8_t mapq2 = proper_pair_mapq(nams_pair[1]);
            output_hits_paf_PE(outstring, nam_read1, record1.name,
                            references,
                            record1.seq.length(), mapq1);
            output_hits_paf_PE(outstring, nam_read2, record2.name,
                            references,
                            record2.seq.length(), mapq2);
        }
    } else {
        Read read1(record1.seq);
        Read read2(record2.seq);
        auto alignment_pairs = align_paired(
            aligner, nams_pair[0], nams_pair[1], read1, read2,
            index_parameters.syncmer.k, references, details,
            effective_params.dropoff_threshold, isize_est,
            effective_params.max_tries,
            effective_params.rescue_sigma_mult
        );

        // -1 marks the typical case that both reads map uniquely and form a
        // proper pair. Then the mapping quality is computed based on the NAMs.
        if (alignment_pairs.size() == 1 && alignment_pairs[0].score == -1) {
            details[0].xp = 5;
            details[1].xp = 5;
            details[0].yj = 1.0f;  // Only one pair on fast path
            details[1].yj = 1.0f;
            Alignment& alignment1 = alignment_pairs[0].alignment1;
            Alignment& alignment2 = alignment_pairs[0].alignment2;
            bool is_proper = is_proper_pair(alignment1, alignment2, isize_est.mu, isize_est.sigma);
            if (
                is_proper
                && isize_est.sample_size < 400
                && alignment1.edit_distance + alignment2.edit_distance < 3
            ) {
                isize_est.update(std::abs(alignment1.ref_start - alignment2.ref_start));
            }

            uint8_t mapq1 = proper_pair_mapq(nams_pair[0]);
            uint8_t mapq2 = proper_pair_mapq(nams_pair[1]);

            details[0].best_alignments = 1;
            details[1].best_alignments = 1;

            // Build SV extra tags if requested (skip in pass 1 — output is suppressed)
            std::string sv_extra1, sv_extra2;
            if (map_param.sv_tags && !collector) {
                sv_extra1 = build_sv_tags(alignment1, alignment2, isize_est.mu, isize_est.sigma, details[0], true);
                sv_extra2 = build_sv_tags(alignment1, alignment2, isize_est.mu, isize_est.sigma, details[1], true);
            }

            if (effective_params.max_supplementary > 0) {
                auto k = index_parameters.syncmer.k;
                auto si1 = build_supplementary_info(sam, aligner, nams_pair[0], read1, alignment1, mapq1, k, references, effective_params.dropoff_threshold, effective_params.max_tries, effective_params.max_supplementary, effective_params.max_supp_overlap, effective_params.min_clip);
                auto si2 = build_supplementary_info(sam, aligner, nams_pair[1], read2, alignment2, mapq2, k, references, effective_params.dropoff_threshold, effective_params.max_tries, effective_params.max_supplementary, effective_params.max_supp_overlap, effective_params.min_clip);

                // Combine SA tags with SV tags
                std::string extra1 = si1.primary_sa_tag;
                std::string extra2 = si2.primary_sa_tag;
                if (!sv_extra1.empty()) {
                    if (!extra1.empty()) extra1 += '\t';
                    extra1 += sv_extra1;
                }
                if (!sv_extra2.empty()) {
                    if (!extra2.empty()) extra2 += '\t';
                    extra2 += sv_extra2;
                }
                // Collect SV evidence for two-pass (pass 1) — use cheap position
                // extraction from supplementary alignments, skip expensive breakpoint
                // computation (microhomology, CI width, formatting)
                if (collector) {
                    for (const auto& s : si1.supps) {
                        collector->add(alignment1.ref_id, alignment1.ref_start + alignment1.length);
                        collector->add(s.ref_id, s.ref_start);
                    }
                    for (const auto& s : si2.supps) {
                        collector->add(alignment2.ref_id, alignment2.ref_start + alignment2.length);
                        collector->add(s.ref_id, s.ref_start);
                    }
                    if (!is_proper && !alignment1.is_unaligned && !alignment2.is_unaligned
                        && isize_est.sample_size >= 100) {
                        collector->add(alignment1.ref_id, alignment1.ref_start);
                        collector->add(alignment2.ref_id, alignment2.ref_start);
                    }
                }

                // XD/XR + breakpoint tags (pass 2 / normal mode only)
                std::vector<BreakpointInfo> bp_infos1, bp_infos2;
                if (map_param.sv_tags && !collector) {
                    if (!si1.supps.empty()) {
                        auto dt = compute_supp_delta_tags(alignment1, si1.supps);
                        if (!dt.empty()) {
                            if (!extra1.empty()) extra1 += '\t';
                            extra1 += dt;
                        }
                        for (const auto& s : si1.supps) {
                            bp_infos1.push_back(compute_breakpoint_info(alignment1, s, references, record1.seq));
                        }
                        if (!bp_infos1.empty()) {
                            // YS on primary: comma-separated SV types
                            if (!extra1.empty()) extra1 += '\t';
                            extra1 += "YS:Z:";
                            for (size_t j = 0; j < bp_infos1.size(); ++j) {
                                if (j > 0) extra1 += ',';
                                extra1 += bp_infos1[j].sv_type;
                            }
                            extra1 += '\t';
                            extra1 += bp_infos1[0].format_all(references);
                        }
                    }
                    if (!si2.supps.empty()) {
                        auto dt = compute_supp_delta_tags(alignment2, si2.supps);
                        if (!dt.empty()) {
                            if (!extra2.empty()) extra2 += '\t';
                            extra2 += dt;
                        }
                        for (const auto& s : si2.supps) {
                            bp_infos2.push_back(compute_breakpoint_info(alignment2, s, references, record2.seq));
                        }
                        if (!bp_infos2.empty()) {
                            // YS on primary: comma-separated SV types
                            if (!extra2.empty()) extra2 += '\t';
                            extra2 += "YS:Z:";
                            for (size_t j = 0; j < bp_infos2.size(); ++j) {
                                if (j > 0) extra2 += ',';
                                extra2 += bp_infos2[j].sv_type;
                            }
                            extra2 += '\t';
                            extra2 += bp_infos2[0].format_all(references);
                        }
                    }
                }

                // XK: soft-clip realignment for short clips
                if (map_param.sv_tags && !collector) {
                    bool lc1 = false, rc1 = false, lc2 = false, rc2 = false;
                    for (const auto& s : si1.supps) {
                        if (s.query_start < alignment1.query_start) lc1 = true;
                        if (s.query_end > alignment1.query_end) rc1 = true;
                    }
                    for (const auto& s : si2.supps) {
                        if (s.query_start < alignment2.query_start) lc2 = true;
                        if (s.query_end > alignment2.query_end) rc2 = true;
                    }
                    if (alignment1.is_revcomp) std::swap(lc1, rc1);
                    if (alignment2.is_revcomp) std::swap(lc2, rc2);
                    const auto& seq1 = alignment1.is_revcomp ? read1.rc : read1.seq;
                    auto xk1 = realign_soft_clips(alignment1, seq1, references, aligner, lc1, rc1);
                    if (!xk1.empty()) {
                        if (!extra1.empty()) extra1 += '\t';
                        extra1 += "XK:Z:" + format_xk_tag(xk1);
                    }
                    const auto& seq2 = alignment2.is_revcomp ? read2.rc : read2.seq;
                    auto xk2 = realign_soft_clips(alignment2, seq2, references, aligner, lc2, rc2);
                    if (!xk2.empty()) {
                        if (!extra2.empty()) extra2 += '\t';
                        extra2 += "XK:Z:" + format_xk_tag(xk2);
                    }
                }

                // Tag hotspot-realigned reads
                if (is_near_hotspot) {
                    if (!extra1.empty()) extra1 += '\t';
                    extra1 += "XH:i:1";
                    if (!extra2.empty()) extra2 += '\t';
                    extra2 += "XH:i:1";
                }

                if (phasing_map) {
                    if (!alignment1.is_unaligned) {
                        auto p1 = phasing_map->assign_haplotype(alignment1, record1.seq, alignment1.cigar);
                        if (p1) { if (!extra1.empty()) extra1 += '\t'; extra1 += format_phase_tags(*p1); }
                    }
                    if (!alignment2.is_unaligned) {
                        auto p2 = phasing_map->assign_haplotype(alignment2, record2.seq, alignment2.cigar);
                        if (p2) { if (!extra2.empty()) extra2 += '\t'; extra2 += format_phase_tags(*p2); }
                    }
                }

                sam.add_pair(alignment1, alignment2, record1, record2, read1.rc, read2.rc, mapq1, mapq2, is_proper, PRIMARY, details, extra1, extra2);

                for (size_t j = 0; j < si1.supps.size(); ++j) {
                    std::string supp_extra = si1.supp_sa_tags[j];
                    if (map_param.sv_tags) {
                        if (!supp_extra.empty()) supp_extra += '\t';
                        supp_extra += compute_ys_tag(si1.supps[j], alignment1);
                        if (j < bp_infos1.size()) {
                            supp_extra += '\t';
                            supp_extra += bp_infos1[j].format_all(references);
                        }
                        if (j < si1.tied_counts.size() && si1.tied_counts[j] > 0) {
                            supp_extra += '\t';
                            supp_extra += "Xt:i:" + std::to_string(si1.tied_counts[j]);
                        }
                    }
                    sam.add_paired_supplementary(si1.supps[j], record1, read1.rc, mapq1, true, alignment2, mapq2, details[0], supp_extra);
                }
                for (size_t j = 0; j < si2.supps.size(); ++j) {
                    std::string supp_extra = si2.supp_sa_tags[j];
                    if (map_param.sv_tags) {
                        if (!supp_extra.empty()) supp_extra += '\t';
                        supp_extra += compute_ys_tag(si2.supps[j], alignment2);
                        if (j < bp_infos2.size()) {
                            supp_extra += '\t';
                            supp_extra += bp_infos2[j].format_all(references);
                        }
                        if (j < si2.tied_counts.size() && si2.tied_counts[j] > 0) {
                            supp_extra += '\t';
                            supp_extra += "Xt:i:" + std::to_string(si2.tied_counts[j]);
                        }
                    }
                    sam.add_paired_supplementary(si2.supps[j], record2, read2.rc, mapq2, false, alignment1, mapq1, details[1], supp_extra);
                }
            } else {
                if (phasing_map) {
                    if (!alignment1.is_unaligned) {
                        auto p1 = phasing_map->assign_haplotype(alignment1, record1.seq, alignment1.cigar);
                        if (p1) { if (!sv_extra1.empty()) sv_extra1 += '\t'; sv_extra1 += format_phase_tags(*p1); }
                    }
                    if (!alignment2.is_unaligned) {
                        auto p2 = phasing_map->assign_haplotype(alignment2, record2.seq, alignment2.cigar);
                        if (p2) { if (!sv_extra2.empty()) sv_extra2 += '\t'; sv_extra2 += format_phase_tags(*p2); }
                    }
                }
                sam.add_pair(alignment1, alignment2, record1, record2, read1.rc, read2.rc, mapq1, mapq2, is_proper, PRIMARY, details, sv_extra1, sv_extra2);
            }
        } else {
            std::sort(alignment_pairs.begin(), alignment_pairs.end(), by_score<ScoredAlignmentPair>);
            deduplicate_scored_pairs(alignment_pairs);

            // If there are multiple top-scoring alignments (all with the same score),
            // pick one randomly and move it to the front.
            size_t i = count_best_alignment_pairs(alignment_pairs);
            details[0].best_alignments = i;
            details[1].best_alignments = i;
            if (i > 1) {
                size_t random_index = std::uniform_int_distribution<>(0, i - 1)(random_engine);
                if (random_index != 0) {
                    std::swap(alignment_pairs[0], alignment_pairs[random_index]);
                }
            }

            // Set XS (best non-primary alignment score per read) and YJ (pair posterior)
            if (alignment_pairs.size() > 1) {
                auto& second = alignment_pairs[1];
                if (!second.alignment1.is_unaligned) details[0].xs = second.alignment1.score;
                if (!second.alignment2.is_unaligned) details[1].xs = second.alignment2.score;
            }
            {
                float posterior = compute_pair_posterior(alignment_pairs);
                details[0].yj = posterior;
                details[1].yj = posterior;
            }

            if (effective_params.max_supplementary > 0 && !alignment_pairs.empty()) {
                // Detect supplementaries for each read, then output
                // primary pair with SA tags, supplementaries, and secondaries
                auto [mapq1, mapq2] = joint_mapq_from_high_scores(alignment_pairs);
                auto& primary1 = alignment_pairs[0].alignment1;
                auto& primary2 = alignment_pairs[0].alignment2;
                bool is_proper = is_proper_pair(primary1, primary2, isize_est.mu, isize_est.sigma);

                // Extract unique per-read alignments from scored pairs
                // (these were already computed by align_paired, no need to re-align)
                std::vector<Alignment> candidates1, candidates2;
                for (const auto& ap : alignment_pairs) {
                    const auto& a1 = ap.alignment1;
                    if (!a1.is_unaligned) {
                        bool found = false;
                        for (const auto& c : candidates1) {
                            if (c.ref_id == a1.ref_id && c.ref_start == a1.ref_start && c.is_revcomp == a1.is_revcomp) { found = true; break; }
                        }
                        if (!found) candidates1.push_back(a1);
                    }
                    const auto& a2 = ap.alignment2;
                    if (!a2.is_unaligned) {
                        bool found = false;
                        for (const auto& c : candidates2) {
                            if (c.ref_id == a2.ref_id && c.ref_start == a2.ref_start && c.is_revcomp == a2.is_revcomp) { found = true; break; }
                        }
                        if (!found) candidates2.push_back(a2);
                    }
                }

                auto si1 = build_supplementary_info_from_candidates(sam, candidates1, primary1, mapq1, effective_params.max_supplementary, effective_params.max_supp_overlap, record1.seq.length(), effective_params.min_clip);
                auto si2 = build_supplementary_info_from_candidates(sam, candidates2, primary2, mapq2, effective_params.max_supplementary, effective_params.max_supp_overlap, record2.seq.length(), effective_params.min_clip);

                // Build SV extra tags if requested (skip in pass 1 — output is suppressed)
                std::string sv_extra1, sv_extra2;
                if (map_param.sv_tags && !collector) {
                    sv_extra1 = build_sv_tags(primary1, primary2, isize_est.mu, isize_est.sigma, details[0], false);
                    sv_extra2 = build_sv_tags(primary1, primary2, isize_est.mu, isize_est.sigma, details[1], false);
                }

                // Combine SA tags with SV tags
                std::string extra1 = si1.primary_sa_tag;
                std::string extra2 = si2.primary_sa_tag;
                if (!sv_extra1.empty()) {
                    if (!extra1.empty()) extra1 += '\t';
                    extra1 += sv_extra1;
                }
                if (!sv_extra2.empty()) {
                    if (!extra2.empty()) extra2 += '\t';
                    extra2 += sv_extra2;
                }
                // Collect SV evidence for two-pass (pass 1) — cheap position extraction
                if (collector) {
                    for (const auto& s : si1.supps) {
                        collector->add(primary1.ref_id, primary1.ref_start + primary1.length);
                        collector->add(s.ref_id, s.ref_start);
                    }
                    for (const auto& s : si2.supps) {
                        collector->add(primary2.ref_id, primary2.ref_start + primary2.length);
                        collector->add(s.ref_id, s.ref_start);
                    }
                    if (!is_proper && !primary1.is_unaligned && !primary2.is_unaligned
                        && isize_est.sample_size >= 100) {
                        collector->add(primary1.ref_id, primary1.ref_start);
                        collector->add(primary2.ref_id, primary2.ref_start);
                    }
                }

                // XD/XR + breakpoint tags (pass 2 / normal mode only)
                std::vector<BreakpointInfo> bp_infos1, bp_infos2;
                if (map_param.sv_tags && !collector) {
                    if (!si1.supps.empty()) {
                        auto dt = compute_supp_delta_tags(primary1, si1.supps);
                        if (!dt.empty()) {
                            if (!extra1.empty()) extra1 += '\t';
                            extra1 += dt;
                        }
                        for (const auto& s : si1.supps) {
                            bp_infos1.push_back(compute_breakpoint_info(primary1, s, references, record1.seq));
                        }
                        if (!bp_infos1.empty()) {
                            // YS on primary: comma-separated SV types
                            if (!extra1.empty()) extra1 += '\t';
                            extra1 += "YS:Z:";
                            for (size_t j = 0; j < bp_infos1.size(); ++j) {
                                if (j > 0) extra1 += ',';
                                extra1 += bp_infos1[j].sv_type;
                            }
                            extra1 += '\t';
                            extra1 += bp_infos1[0].format_all(references);
                        }
                    }
                    if (!si2.supps.empty()) {
                        auto dt = compute_supp_delta_tags(primary2, si2.supps);
                        if (!dt.empty()) {
                            if (!extra2.empty()) extra2 += '\t';
                            extra2 += dt;
                        }
                        for (const auto& s : si2.supps) {
                            bp_infos2.push_back(compute_breakpoint_info(primary2, s, references, record2.seq));
                        }
                        if (!bp_infos2.empty()) {
                            // YS on primary: comma-separated SV types
                            if (!extra2.empty()) extra2 += '\t';
                            extra2 += "YS:Z:";
                            for (size_t j = 0; j < bp_infos2.size(); ++j) {
                                if (j > 0) extra2 += ',';
                                extra2 += bp_infos2[j].sv_type;
                            }
                            extra2 += '\t';
                            extra2 += bp_infos2[0].format_all(references);
                        }
                    }
                }

                // XK: soft-clip realignment for short clips
                if (map_param.sv_tags && !collector) {
                    bool lc1 = false, rc1 = false, lc2 = false, rc2 = false;
                    for (const auto& s : si1.supps) {
                        if (s.query_start < primary1.query_start) lc1 = true;
                        if (s.query_end > primary1.query_end) rc1 = true;
                    }
                    for (const auto& s : si2.supps) {
                        if (s.query_start < primary2.query_start) lc2 = true;
                        if (s.query_end > primary2.query_end) rc2 = true;
                    }
                    if (primary1.is_revcomp) std::swap(lc1, rc1);
                    if (primary2.is_revcomp) std::swap(lc2, rc2);
                    const auto& seq1 = primary1.is_revcomp ? read1.rc : read1.seq;
                    auto xk1 = realign_soft_clips(primary1, seq1, references, aligner, lc1, rc1);
                    if (!xk1.empty()) {
                        if (!extra1.empty()) extra1 += '\t';
                        extra1 += "XK:Z:" + format_xk_tag(xk1);
                    }
                    const auto& seq2 = primary2.is_revcomp ? read2.rc : read2.seq;
                    auto xk2 = realign_soft_clips(primary2, seq2, references, aligner, lc2, rc2);
                    if (!xk2.empty()) {
                        if (!extra2.empty()) extra2 += '\t';
                        extra2 += "XK:Z:" + format_xk_tag(xk2);
                    }
                }

                // Tag hotspot-realigned reads
                if (is_near_hotspot) {
                    if (!extra1.empty()) extra1 += '\t';
                    extra1 += "XH:i:1";
                    if (!extra2.empty()) extra2 += '\t';
                    extra2 += "XH:i:1";
                }

                // HP/PS: haplotype phasing
                if (phasing_map) {
                    if (!primary1.is_unaligned) {
                        auto p1 = phasing_map->assign_haplotype(primary1, record1.seq, primary1.cigar);
                        if (p1) { if (!extra1.empty()) extra1 += '\t'; extra1 += format_phase_tags(*p1); }
                    }
                    if (!primary2.is_unaligned) {
                        auto p2 = phasing_map->assign_haplotype(primary2, record2.seq, primary2.cigar);
                        if (p2) { if (!extra2.empty()) extra2 += '\t'; extra2 += format_phase_tags(*p2); }
                    }
                }

                // Output primary pair with SA + SV tags
                sam.add_pair(primary1, primary2, record1, record2, read1.rc, read2.rc, mapq1, mapq2, is_proper, PRIMARY, details, extra1, extra2);

                // Output supplementary records
                for (size_t j = 0; j < si1.supps.size(); ++j) {
                    std::string supp_extra = si1.supp_sa_tags[j];
                    if (map_param.sv_tags) {
                        if (!supp_extra.empty()) supp_extra += '\t';
                        supp_extra += compute_ys_tag(si1.supps[j], primary1);
                        if (j < bp_infos1.size()) {
                            supp_extra += '\t';
                            supp_extra += bp_infos1[j].format_all(references);
                        }
                        if (j < si1.tied_counts.size() && si1.tied_counts[j] > 0) {
                            supp_extra += '\t';
                            supp_extra += "Xt:i:" + std::to_string(si1.tied_counts[j]);
                        }
                    }
                    sam.add_paired_supplementary(si1.supps[j], record1, read1.rc, mapq1, true, primary2, mapq2, details[0], supp_extra);
                }
                for (size_t j = 0; j < si2.supps.size(); ++j) {
                    std::string supp_extra = si2.supp_sa_tags[j];
                    if (map_param.sv_tags) {
                        if (!supp_extra.empty()) supp_extra += '\t';
                        supp_extra += compute_ys_tag(si2.supps[j], primary2);
                        if (j < bp_infos2.size()) {
                            supp_extra += '\t';
                            supp_extra += bp_infos2[j].format_all(references);
                        }
                        if (j < si2.tied_counts.size() && si2.tied_counts[j] > 0) {
                            supp_extra += '\t';
                            supp_extra += "Xt:i:" + std::to_string(si2.tied_counts[j]);
                        }
                    }
                    sam.add_paired_supplementary(si2.supps[j], record2, read2.rc, mapq2, false, primary1, mapq1, details[1], supp_extra);
                }

                // Output secondary pairs
                double secondary_dropoff = 2 * aligner.parameters.mismatch + aligner.parameters.gap_open;
                auto max_out = std::min(alignment_pairs.size(), static_cast<size_t>(map_param.max_secondary) + 1);
                auto s_max = alignment_pairs[0].score;
                for (size_t j = 1; j < max_out; ++j) {
                    auto& pair = alignment_pairs[j];
                    if (s_max - pair.score >= secondary_dropoff) break;
                    bool sec_proper = is_proper_pair(pair.alignment1, pair.alignment2, isize_est.mu, isize_est.sigma);
                    sam.add_pair(pair.alignment1, pair.alignment2, record1, record2, read1.rc, read2.rc, 0, 0, sec_proper, SECONDARY_ALN, details);
                }
            } else {
                double secondary_dropoff = 2 * aligner.parameters.mismatch + aligner.parameters.gap_open;
                output_aligned_pairs(
                    alignment_pairs,
                    sam,
                    map_param.max_secondary,
                    secondary_dropoff,
                    record1,
                    record2,
                    read1,
                    read2,
                    isize_est.mu,
                    isize_est.sigma,
                    details,
                    map_param.sv_tags,
                    phasing_map
                );
            }
        }
    }
    statistics.tot_extend += extend_timer.duration();
    statistics += details[0];
    statistics += details[1];
}

void align_or_map_single(
    const KSeq &record,
    Sam& sam,
    std::string &outstring,
    AlignmentStatistics &statistics,
    const Aligner &aligner,
    const Chainer& chainer,
    const MappingParameters &map_param,
    const IndexParameters& index_parameters,
    const References& references,
    const StrobemerIndex& index,
    std::minstd_rand& random_engine,
    std::vector<double> &abundances,
    SvEvidenceCollector* collector,
    const HotspotMap* hotspot_map,
    const MappingParameters* relaxed_params,
    const PhasingMap* phasing_map
) {
    Details details;
    std::vector<Nam> nams;

    logger.trace() << "\nQuery: " << record.name << '\n';
    nams = get_nams_or_chains(record, index, chainer, statistics, details, map_param, index_parameters, random_engine);
    if (!nams.empty()) {
        details.s1 = std::lroundf(nams[0].score);
        if (nams.size() > 1) details.s2 = std::lroundf(nams[1].score);
        details.cm = nams[0].n_matches;
        details.xc = std::lroundf(nams[0].score);
    }

    Timer extend_timer;
    size_t n_best = 0;
    switch (map_param.output_format) {
        case OutputFormat::Abundance: {
            if (!nams.empty()){
                for (auto &t : nams){
                    if (t.score == nams[0].score){
                        ++n_best;
                    }else{
                        break;
                    }
                }

                for (auto &nam: nams) {
                    if (nam.ref_start < 0) {
                        continue;
                    }
                    if (nam.score != nams[0].score){
                        break;
                    }
                    abundances[nam.ref_id] += float(record.seq.length()) / float(n_best);
                }
            }
        }
        break;
        case OutputFormat::PAF: {
            int mapq = proper_pair_mapq(nams);
            output_hits_paf(outstring, nams, record.name, references, record.seq.length(), mapq);
            break;
        }
        case OutputFormat::SAM: {
            // Two-pass: check if read's top NAM is near an SV hotspot
            bool se_near_hotspot = false;
            if (hotspot_map && !nams.empty()) {
                se_near_hotspot = hotspot_map->near_hotspot(nams[0].ref_id, nams[0].ref_start);
            }
            const auto& se_params = (se_near_hotspot && relaxed_params) ? *relaxed_params : map_param;
            align_single(
                aligner, sam, nams, record, index_parameters.syncmer.k,
                references, details, se_params.dropoff_threshold, se_params.max_tries,
                se_params.max_secondary, se_params.max_supplementary,
                se_params.max_supp_overlap, random_engine, se_params.sv_tags,
                collector, se_near_hotspot, phasing_map, se_params.min_clip
            );
            break;
        }
    }
    statistics.tot_extend += extend_timer.duration();
    statistics += details;
}

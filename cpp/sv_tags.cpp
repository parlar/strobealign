#include "sv_tags.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

std::string BreakpointInfo::format_all(const References& refs) const {
    std::string result;
    // YB: breakpoint intervals
    int left_lo = std::max(0, left_bp - mh_bwd);
    int left_hi = left_bp + mh_fwd;
    int right_lo = std::max(0, right_bp - mh_bwd);
    int right_hi = right_bp + mh_fwd;
    result += "YB:Z:";
    result += refs.names[left_ref_id];
    result += ':';
    result += std::to_string(left_lo);
    result += '-';
    result += std::to_string(left_hi);
    result += '|';
    result += refs.names[right_ref_id];
    result += ':';
    result += std::to_string(right_lo);
    result += '-';
    result += std::to_string(right_hi);
    // YC: CI widths
    result += "\tYC:Z:";
    result += std::to_string(left_ci);
    result += ',';
    result += std::to_string(right_ci);
    // YM: microhomology
    result += "\tYM:i:";
    result += std::to_string(mh_fwd + mh_bwd);
    // YI: inserted sequence (only if present)
    if (!inserted_seq.empty()) {
        result += "\tYI:Z:";
        result += inserted_seq;
    }
    // YA: alternative placements (only if microhomology > 0 and <= 20)
    int mh_total = mh_fwd + mh_bwd;
    if (mh_total > 0 && mh_total <= 20) {
        result += "\tYA:Z:";
        bool first = true;
        for (int i = -mh_bwd; i <= mh_fwd; ++i) {
            if (!first) result += ';';
            result += std::to_string(left_bp + i);
            result += ',';
            result += std::to_string(right_bp + i);
            first = false;
        }
    }
    return result;
}

std::string compute_yt_tag(
    const Alignment& a1, const Alignment& a2, float mu, float sigma
) {
    if (a1.is_unaligned && a2.is_unaligned) {
        return "YT:Z:UU";
    }
    if (a1.is_unaligned || a2.is_unaligned) {
        return "";
    }
    if (a1.ref_id != a2.ref_id) {
        return "YT:Z:IC";
    }
    // Same chromosome, same strand → deviant orientation
    if (a1.is_revcomp == a2.is_revcomp) {
        return "YT:Z:DO";
    }
    // Opposite strands — use is_proper_pair for CP to match PROPER_PAIR flag
    if (is_proper_pair(a1, a2, mu, sigma)) {
        return "YT:Z:CP";
    }
    return "YT:Z:DI";
}

/*
 * Compute YD insert-size Z-score: |insert_size - mu| / sigma
 * Only meaningful when both reads map to the same chromosome.
 */
std::string compute_yd_tag(
    const Alignment& a1, const Alignment& a2, float mu, float sigma
) {
    if (a1.is_unaligned || a2.is_unaligned) return "";
    if (a1.ref_id != a2.ref_id) return "";
    if (sigma <= 0) return "";
    float insert = std::abs(a2.ref_start - a1.ref_start);
    float zscore = std::fabs(insert - mu) / sigma;
    char buf[32];
    std::snprintf(buf, sizeof(buf), "YD:f:%.2f", zscore);
    return buf;
}

/*
 * Compute Xr rescue-status tag:
 *   0 - full search (normal alignment, not fast path and no rescue)
 *   1 - fast path (score == -1 in alignment pairs)
 *   2 - mate was rescued
 */
std::string compute_xr_tag(const Details& details, bool is_fast_path) {
    if (details.mate_rescue > 0) {
        return "Xr:i:2";
    }
    if (is_fast_path) {
        return "Xr:i:1";
    }
    return "Xr:i:0";
}

/*
 * Classify SV type for supplementary vs primary:
 *   TRA - different chromosomes (translocation)
 *   INV - same chromosome, different strands (inversion)
 *   DUP - same chromosome, same strand, reference ranges overlap (duplication)
 *   DEL - same chromosome, same strand, no reference overlap (deletion)
 */
std::string classify_sv_type(const Alignment& supp, const Alignment& primary) {
    if (supp.ref_id != primary.ref_id) {
        return "TRA";
    }
    if (supp.is_revcomp != primary.is_revcomp) {
        return "INV";
    }
    int supp_end = supp.ref_start + supp.length;
    int prim_end = primary.ref_start + primary.length;
    bool overlap = (supp.ref_start < prim_end && primary.ref_start < supp_end);
    return overlap ? "DUP" : "DEL";
}

std::string compute_ys_tag(const Alignment& supp, const Alignment& primary) {
    return "YS:Z:" + classify_sv_type(supp, primary);
}

int facing_ref_pos(const Alignment& aln, bool covers_left_of_read) {
    if (covers_left_of_read) {
        return aln.is_revcomp ? aln.ref_start : (aln.ref_start + aln.length);
    } else {
        return aln.is_revcomp ? (aln.ref_start + aln.length) : aln.ref_start;
    }
}

/*
 * Compute microhomology by bidirectional reference scan at breakpoint junction.
 * Returns {forward_matches, backward_matches}.
 */
std::pair<int, int> compute_microhomology(
    const References& references,
    int left_ref_id, int left_bp,
    int right_ref_id, int right_bp
) {
    const auto& left_seq = references.sequences[left_ref_id];
    const auto& right_seq = references.sequences[right_ref_id];

    // Forward scan
    int fwd = 0;
    int max_fwd = std::min({
        static_cast<int>(left_seq.size()) - left_bp,
        static_cast<int>(right_seq.size()) - right_bp,
        50
    });
    for (int i = 0; i < max_fwd; ++i) {
        char l = std::toupper(left_seq[left_bp + i]);
        char r = std::toupper(right_seq[right_bp + i]);
        if (l == 'N' || r == 'N' || l != r) break;
        fwd++;
    }

    // Backward scan
    int bwd = 0;
    int max_bwd = std::min({left_bp, right_bp, 50});
    for (int i = 0; i < max_bwd; ++i) {
        char l = std::toupper(left_seq[left_bp - 1 - i]);
        char r = std::toupper(right_seq[right_bp - 1 - i]);
        if (l == 'N' || r == 'N' || l != r) break;
        bwd++;
    }
    return {fwd, bwd};
}

/*
 * Detect inserted sequence at the breakpoint junction.
 * Uses forward-strand query coordinates to find uncovered read bases.
 */
std::string compute_inserted_sequence(
    const Alignment& primary,
    const Alignment& supp,
    const std::string& read_seq
) {
    if (primary.query_end <= supp.query_start) {
        int gap_start = primary.query_end;
        int gap_end = supp.query_start;
        if (gap_end > gap_start && gap_end <= static_cast<int>(read_seq.size())) {
            return read_seq.substr(gap_start, gap_end - gap_start);
        }
    } else if (supp.query_end <= primary.query_start) {
        int gap_start = supp.query_end;
        int gap_end = primary.query_start;
        if (gap_end > gap_start && gap_end <= static_cast<int>(read_seq.size())) {
            return read_seq.substr(gap_start, gap_end - gap_start);
        }
    }
    return "";
}

/*
 * Compute confidence interval width for a breakpoint position.
 * CI = microhomology_total + round(local_repeat_fraction * 10)
 */
int compute_ci_width(
    const References& references,
    int ref_id, int bp_pos,
    int microhomology_total
) {
    const auto& seq = references.sequences[ref_id];
    int window = 50;
    int start = std::max(0, bp_pos - window);
    int end = std::min(static_cast<int>(seq.size()), bp_pos + window);
    int total_bases = end - start;
    if (total_bases <= 0) return microhomology_total;

    // Mark repeat positions to avoid double-counting between homopolymer and
    // dinucleotide repeat detectors
    std::vector<bool> is_repeat(total_bases, false);

    // Detect homopolymer runs (>=3bp)
    for (int i = start; i < end; ) {
        char c = std::toupper(seq[i]);
        int run = 1;
        while (i + run < end && std::toupper(seq[i + run]) == c) run++;
        if (run >= 3) {
            for (int j = 0; j < run; ++j) is_repeat[i - start + j] = true;
        }
        i += run;
    }

    // Detect dinucleotide repeats (>=6bp = 3 units)
    for (int i = start; i < end - 5; ) {
        char c1 = std::toupper(seq[i]);
        char c2 = std::toupper(seq[i + 1]);
        if (c1 == c2) { i++; continue; }
        int run = 2;
        while (i + run + 1 < end &&
               std::toupper(seq[i + run]) == c1 &&
               std::toupper(seq[i + run + 1]) == c2) {
            run += 2;
        }
        if (run >= 6) {
            for (int j = 0; j < run; ++j) is_repeat[i - start + j] = true;
        }
        i += std::max(run, 1);
    }

    int repeat_bases = 0;
    for (bool r : is_repeat) repeat_bases += r;
    float repeat_frac = static_cast<float>(repeat_bases) / total_bases;
    int repeat_bonus = std::lround(repeat_frac * 10.0f);
    return microhomology_total + repeat_bonus;
}

BreakpointInfo compute_breakpoint_info(
    const Alignment& primary,
    const Alignment& supp,
    const References& references,
    const std::string& read_seq
) {
    BreakpointInfo bp;
    bp.sv_type = classify_sv_type(supp, primary);

    if (bp.sv_type == "DEL") {
        if (primary.ref_start <= supp.ref_start) {
            bp.left_bp = primary.ref_start + primary.length;
            bp.right_bp = supp.ref_start;
        } else {
            bp.left_bp = supp.ref_start + supp.length;
            bp.right_bp = primary.ref_start;
        }
        bp.left_ref_id = bp.right_ref_id = primary.ref_id;
    } else if (bp.sv_type == "DUP") {
        if (primary.ref_start <= supp.ref_start) {
            bp.left_bp = supp.ref_start;
            bp.right_bp = primary.ref_start + primary.length;
        } else {
            bp.left_bp = primary.ref_start;
            bp.right_bp = supp.ref_start + supp.length;
        }
        bp.left_ref_id = bp.right_ref_id = primary.ref_id;
    } else {
        bool primary_covers_left = (primary.query_start <= supp.query_start);
        const auto& left_aln = primary_covers_left ? primary : supp;
        const auto& right_aln = primary_covers_left ? supp : primary;

        bp.left_bp = facing_ref_pos(left_aln, true);
        bp.right_bp = facing_ref_pos(right_aln, false);
        bp.left_ref_id = left_aln.ref_id;
        bp.right_ref_id = right_aln.ref_id;
    }

    bp.left_bp = std::max(0, std::min(bp.left_bp, static_cast<int>(references.sequences[bp.left_ref_id].size())));
    bp.right_bp = std::max(0, std::min(bp.right_bp, static_cast<int>(references.sequences[bp.right_ref_id].size())));

    if (bp.left_ref_id == bp.right_ref_id && bp.left_bp == bp.right_bp) {
        bp.mh_fwd = 0;
        bp.mh_bwd = 0;
    } else {
        auto [fwd, bwd] = compute_microhomology(references, bp.left_ref_id, bp.left_bp, bp.right_ref_id, bp.right_bp);
        bp.mh_fwd = fwd;
        bp.mh_bwd = bwd;
    }

    bp.inserted_seq = compute_inserted_sequence(primary, supp, read_seq);

    int mh_total = bp.mh_fwd + bp.mh_bwd;
    bp.left_ci = compute_ci_width(references, bp.left_ref_id, bp.left_bp, mh_total);
    bp.right_ci = compute_ci_width(references, bp.right_ref_id, bp.right_bp, mh_total);

    return bp;
}

float compute_split_entropy(const std::vector<Alignment>& alignments) {
    if (alignments.size() <= 1) return 0.0f;
    constexpr double T = 10.0;
    int max_score = alignments[0].score;
    for (const auto& a : alignments) {
        if (a.score > max_score) max_score = a.score;
    }
    double Z = 0.0;
    std::vector<double> weights;
    weights.reserve(alignments.size());
    for (const auto& a : alignments) {
        double w = std::exp((a.score - max_score) / T);
        weights.push_back(w);
        Z += w;
    }
    double H = 0.0;
    for (double w : weights) {
        double p = w / Z;
        if (p > 0.0) {
            H -= p * std::log2(p);
        }
    }
    return static_cast<float>(H);
}

float compute_max_artifact_prob(
    const Alignment& primary,
    const std::vector<Alignment>& supplementaries,
    int read_length
) {
    if (supplementaries.empty() || primary.score <= 0 || read_length <= 0) return 0.0f;
    float max_prob = 0.0f;
    for (const auto& s : supplementaries) {
        float len_frac = static_cast<float>(s.query_end - s.query_start) / read_length;
        float score_frac = static_cast<float>(s.score) / primary.score;
        len_frac = std::min(1.0f, std::max(0.0f, len_frac));
        score_frac = std::min(1.0f, std::max(0.0f, score_frac));
        float artifact = 1.0f - (len_frac * score_frac);
        if (artifact > max_prob) max_prob = artifact;
    }
    return max_prob;
}

std::string compute_supp_delta_tags(const Alignment& primary, const std::vector<Alignment>& supplementaries) {
    if (supplementaries.empty() || primary.score <= 0) return "";
    int best_supp_score = supplementaries[0].score;
    int delta = primary.score - best_supp_score;
    float ratio = static_cast<float>(best_supp_score) / primary.score;
    char buf[64];
    std::snprintf(buf, sizeof(buf), "XD:i:%d\tXR:f:%.3f", delta, ratio);
    return buf;
}

std::string build_sv_tags(
    const Alignment& a1, const Alignment& a2,
    float mu, float sigma,
    const Details& details, bool is_fast_path
) {
    std::string tags;
    auto yt = compute_yt_tag(a1, a2, mu, sigma);
    if (!yt.empty()) {
        tags = yt;
    }
    auto yd = compute_yd_tag(a1, a2, mu, sigma);
    if (!yd.empty()) {
        if (!tags.empty()) tags += '\t';
        tags += yd;
    }
    auto xr = compute_xr_tag(details, is_fast_path);
    if (!tags.empty()) tags += '\t';
    tags += xr;
    return tags;
}

/*
 * Compute XU posteriors: probability of REF / SV / RPT classification.
 *
 * Inputs (all from existing SV tags on the SE primary record):
 *   score_ratio  - XR: best_supp_score / primary_score (0-1, higher = supp is strong)
 *   artifact_prob - XF: max artifact probability (0-1, higher = likely artifact)
 *   entropy       - XE: split chain entropy (0+, higher = more ambiguous)
 *   mappability   - XW: local mappability (0-1, lower = more repetitive)
 *   tied_count    - Xt: number of tied supplementary candidates (0+)
 *
 * Uses log-odds + softmax to produce normalized posteriors.
 */
XuPosteriors compute_xu_posteriors(
    float score_ratio,
    float artifact_prob,
    float entropy,
    float mappability,
    int tied_count
) {
    // Clamp inputs to valid ranges
    score_ratio = std::max(0.0f, std::min(1.0f, score_ratio));
    artifact_prob = std::max(0.0f, std::min(1.0f, artifact_prob));
    entropy = std::max(0.0f, entropy);
    mappability = std::max(0.0f, std::min(1.0f, mappability));

    // Log-odds for each class (higher = more likely)
    // SV: strong supplementary, low artifact, low entropy, good mappability
    float lo_sv = 2.0f * score_ratio
                - 3.0f * artifact_prob
                - 1.0f * std::min(entropy, 3.0f)
                + 1.0f * mappability;

    // RPT: low mappability, high entropy, many tied candidates
    float lo_rpt = -2.0f * mappability
                 + 1.5f * std::min(entropy, 3.0f)
                 + 0.5f * std::min(static_cast<float>(tied_count), 10.0f);

    // REF: high artifact probability, weak supplementary
    float lo_ref = 2.0f * artifact_prob
                 - 2.0f * score_ratio
                 + 0.5f * std::min(entropy, 3.0f);

    // Softmax
    float max_lo = std::max({lo_sv, lo_rpt, lo_ref});
    float e_sv  = std::exp(lo_sv - max_lo);
    float e_rpt = std::exp(lo_rpt - max_lo);
    float e_ref = std::exp(lo_ref - max_lo);
    float total = e_sv + e_rpt + e_ref;

    return XuPosteriors{e_ref / total, e_sv / total, e_rpt / total};
}

std::string format_xu_tag(const XuPosteriors& xu) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "XU:Z:REF:%.2f,SV:%.2f,RPT:%.2f", xu.ref, xu.sv, xu.rpt);
    return buf;
}

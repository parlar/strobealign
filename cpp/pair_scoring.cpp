#include "pair_scoring.hpp"

#include <algorithm>
#include <cmath>
#include "aln.hpp"
#include "robin_hood.h"

std::ostream& operator<<(std::ostream& os, const NamPair& nam_pair) {
    os << "NamPair(score=" << nam_pair.score << ", nam1=" << nam_pair.nam1 << ", nam2=" << nam_pair.nam2 << ")";
    return os;
}

/*
 * Compute Boltzmann-weighted posterior probability for the best pair.
 * P(best) = exp(score_best / T) / sum_i(exp(score_i / T)), T = 10.
 */
float compute_pair_posterior(const std::vector<ScoredAlignmentPair>& pairs) {
    if (pairs.empty()) return 0.0f;
    if (pairs.size() == 1) return 1.0f;
    constexpr double T = 10.0;
    double best_score = pairs[0].score;
    double sum = 0.0;
    for (const auto& p : pairs) {
        sum += std::exp((p.score - best_score) / T);
    }
    return static_cast<float>(1.0 / sum);
}

/*
 * Return mapping quality for a read mapped in a proper pair
 */
uint8_t proper_pair_mapq(const std::vector<Nam> &nams) {
    if (nams.size() <= 1) {
        return 60;
    }
    const float s1 = nams[0].score;
    const float s2 = nams[1].score;
    // from minimap2: MAPQ = 40(1−s2/s1) ·min{1,|M|/10} · log s1
    const float min_matches = std::min(nams[0].n_matches / 10.0, 1.0);
    const int uncapped_mapq = 40 * (1 - s2 / s1) * min_matches * log(s1);
    return std::min(uncapped_mapq, 60);
}

/* Compute paired-end mapping score given best alignments (sorted by score) */
std::pair<int, int> joint_mapq_from_high_scores(const std::vector<ScoredAlignmentPair>& pairs) {
    if (pairs.size() <= 1) {
        return std::make_pair(60, 60);
    }
    auto score1 = pairs[0].score;
    auto score2 = pairs[1].score;
    if (score1 == score2) {
        return std::make_pair(0, 0);
    }
    int mapq;
    const int diff = score1 - score2; // (1.0 - (S1 - S2) / S1);
//  float log10_p = diff > 6 ? -6.0 : -diff; // Corresponds to: p_error= 0.1^diff // change in sw score times rough illumina error rate. This is highly heauristic, but so seem most computations of mapq scores
    if (score1 > 0 && score2 > 0) {
        mapq = std::min(60, diff);
//            mapq1 = -10 * log10_p < 60 ? -10 * log10_p : 60;
    } else if (score1 > 0 && score2 <= 0) {
        mapq = 60;
    } else { // both negative SW one is better
        mapq = 1;
    }
    return std::make_pair(mapq, mapq);
}

float normal_pdf(float x, float mu, float sigma)
{
    static const float inv_sqrt_2pi = 0.3989422804014327;
    const float a = (x - mu) / sigma;

    return inv_sqrt_2pi / sigma * std::exp(-0.5f * a * a);
}

std::vector<ScoredAlignmentPair> get_best_scoring_pairs(
    const std::vector<Alignment>& alignments1,
    const std::vector<Alignment>& alignments2,
    float mu,
    float sigma
) {
    // Precompute: log(normal_pdf(x,mu,sigma)) = -0.5*((x-mu)/sigma)^2 - log(sigma) - 0.5*log(2*pi)
    const double inv_sigma = 1.0 / sigma;
    const double log_norm_const = -std::log(sigma) - 0.5 * std::log(2.0 * M_PI);
    const float dist_threshold = mu + 4 * sigma;

    std::vector<ScoredAlignmentPair> pairs;
    for (auto &a1 : alignments1) {
        for (auto &a2 : alignments2) {
            float dist = std::abs(a1.ref_start - a2.ref_start);
            double score = a1.score + a2.score;
            if ((a1.is_revcomp ^ a2.is_revcomp) && (dist < dist_threshold)) {
                double a = (dist - mu) * inv_sigma;
                score += -0.5 * a * a + log_norm_const;
            }
            else { // individual score
                // 10 corresponds to a value of log(normal_pdf(dist, mu, sigma)) of more than 4 stddevs away
                score -= 10;
            }
            pairs.push_back(ScoredAlignmentPair{score, a1, a2});
        }
    }

    return pairs;
}

bool is_proper_nam_pair(const Nam& nam1, const Nam& nam2, float mu, float sigma) {
    if (nam1.ref_id != nam2.ref_id || nam1.is_revcomp == nam2.is_revcomp) {
        return false;
    }
    int r1_ref_start = nam1.projected_ref_start();
    int r2_ref_start = nam2.projected_ref_start();

    // r1 ---> <---- r2
    bool r1_r2 = nam2.is_revcomp && (r1_ref_start <= r2_ref_start) && (r2_ref_start - r1_ref_start < mu + 10*sigma);

     // r2 ---> <---- r1
    bool r2_r1 = nam1.is_revcomp && (r2_ref_start <= r1_ref_start) && (r1_ref_start - r2_ref_start < mu + 10*sigma);

    return r1_r2 || r2_r1;
}

/*
 * Find high-scoring NAMs and NAM pairs. Proper pairs are preferred, but also
 * high-scoring NAMs that could not be paired up are returned (these get a
 * "dummy" NAM as partner in the returned vector).
 */
std::vector<NamPair> get_best_scoring_nam_pairs(
    const std::vector<Nam> &nams1,
    const std::vector<Nam> &nams2,
    float mu,
    float sigma
) {
    std::vector<NamPair> nam_pairs;
    if (nams1.empty() && nams2.empty()) {
        return nam_pairs;
    }

    // Find NAM pairs that appear to be proper pairs
    robin_hood::unordered_set<int> added_n1;
    robin_hood::unordered_set<int> added_n2;
    int best_joint_hits = 0;

    constexpr size_t MAX_NAMS = 1000;
    for (size_t i1 = 0; i1 < std::min(nams1.size(), MAX_NAMS); ++i1) {
        const Nam& nam1 = nams1[i1];
        for (size_t i2 = 0; i2 < std::min(nams2.size(), MAX_NAMS); ++i2) {
            const Nam& nam2 = nams2[i2];
            int joint_hits = nam1.n_matches + nam2.n_matches;
            if (joint_hits < best_joint_hits / 2) {
                break;
            }
            if (is_proper_nam_pair(nam1, nam2, mu, sigma)) {
                nam_pairs.push_back(NamPair{nam1.score + nam2.score, nam1, nam2});
                added_n1.insert(nam1.nam_id);
                added_n2.insert(nam2.nam_id);
                best_joint_hits = std::max(joint_hits, best_joint_hits);
            }
        }
    }

    // Find high-scoring R1 NAMs that are not part of a proper pair
    Nam dummy_nam;
    dummy_nam.ref_start = -1;
    if (!nams1.empty()) {
        int best_joint_hits1 = best_joint_hits > 0 ? best_joint_hits : nams1[0].n_matches;
        for (auto &nam1 : nams1) {
            if (nam1.n_matches < best_joint_hits1 / 2) {
                break;
            }
            if (added_n1.find(nam1.nam_id) != added_n1.end()) {
                continue;
            }
//            int n1_penalty = std::abs(nam1.query_span() - nam1.ref_span());
            nam_pairs.push_back(NamPair{nam1.score, nam1, dummy_nam});
        }
    }

    // Find high-scoring R2 NAMs that are not part of a proper pair
    if (!nams2.empty()) {
        int best_joint_hits2 = best_joint_hits > 0 ? best_joint_hits : nams2[0].n_matches;
        for (auto &nam2 : nams2) {
            if (nam2.n_matches < best_joint_hits2 / 2) {
                break;
            }
            if (added_n2.find(nam2.nam_id) != added_n2.end()){
                continue;
            }
//            int n2_penalty = std::abs(nam2.query_span() - nam2.ref_span());
            nam_pairs.push_back(NamPair{nam2.score, dummy_nam, nam2});
        }
    }

    std::sort(
        nam_pairs.begin(),
        nam_pairs.end(),
        [](const NamPair& a, const NamPair& b) -> bool { return a.score > b.score; }
    ); // Sort by highest score first

    return nam_pairs;
}

/*
 * Align a read to the reference given the mapping location of its mate.
 */
Alignment rescue_align(
    const Aligner& aligner,
    const Nam &mate_nam,
    const References& references,
    const Read& read,
    float mu,
    float sigma,
    int k,
    float sigma_mult
) {
    Alignment alignment;
    int a, b;
    auto read_len = read.size();

    const std::string& r_tmp = mate_nam.is_revcomp ? read.seq : read.rc;
    if (mate_nam.is_revcomp) {
        a = mate_nam.projected_ref_start() - (mu+sigma_mult*sigma);
        b = mate_nam.projected_ref_start() + read_len/2; // at most half read overlap
    } else {
        a = mate_nam.ref_end + (read_len - mate_nam.query_end) - read_len/2; // at most half read overlap
        b = mate_nam.ref_end + (read_len - mate_nam.query_end) + (mu+sigma_mult*sigma);
    }

    auto ref_len = static_cast<int>(references.lengths[mate_nam.ref_id]);
    auto ref_start = std::max(0, std::min(a, ref_len));
    auto ref_end = std::min(ref_len, std::max(0, b));

    if (ref_end < ref_start + k) {
        alignment.cigar = Cigar();
        alignment.edit_distance = read_len;
        alignment.score = 0;
        alignment.ref_start =  0;
        alignment.is_revcomp = !mate_nam.is_revcomp;
        alignment.ref_id = mate_nam.ref_id;
        alignment.is_unaligned = true;
        return alignment;
    }
    const auto& ref_full = references.sequences[mate_nam.ref_id];
    std::string_view ref_segm(ref_full.data() + ref_start, ref_end - ref_start);

    if (!has_shared_substring(r_tmp, ref_segm, k)) {
        alignment.cigar = Cigar();
        alignment.edit_distance = read_len;
        alignment.score = 0;
        alignment.ref_start =  0;
        alignment.is_revcomp = !mate_nam.is_revcomp;
        alignment.ref_id = mate_nam.ref_id;
        alignment.is_unaligned = true;
        return alignment;
    }
    auto opt_info = aligner.align(r_tmp, ref_segm);
    if (opt_info) {
        auto info = opt_info.value();
        alignment.cigar = info.cigar;
        alignment.edit_distance = info.edit_distance;
        alignment.score = info.sw_score;
        alignment.ref_start = ref_start + info.ref_start;
        alignment.is_revcomp = !mate_nam.is_revcomp;
        alignment.ref_id = mate_nam.ref_id;
        alignment.is_unaligned = info.cigar.empty();
        alignment.length = info.ref_span();

        // Convert query coordinates to forward-strand
        if (alignment.is_revcomp) {
            int read_len = static_cast<int>(read.size());
            alignment.query_start = read_len - info.query_end;
            alignment.query_end = read_len - info.query_start;
        } else {
            alignment.query_start = info.query_start;
            alignment.query_end = info.query_end;
        }
    } else {
        alignment.is_unaligned = true;
        alignment.edit_distance = 100000;
        alignment.ref_start = 0;
        alignment.score = -100000;
    }
    return alignment;
}

/*
 * Remove consecutive identical alignment pairs and leave only the first.
 */
void deduplicate_scored_pairs(std::vector<ScoredAlignmentPair>& pairs) {
    if (pairs.size() < 2) {
        return;
    }
    int prev_ref_start1 = pairs[0].alignment1.ref_start;
    int prev_ref_start2 = pairs[0].alignment2.ref_start;
    int prev_ref_id1 = pairs[0].alignment1.ref_id;
    int prev_ref_id2 = pairs[0].alignment2.ref_id;
    size_t j = 1;
    for (size_t i = 1; i < pairs.size(); i++) {
        int ref_start1 = pairs[i].alignment1.ref_start;
        int ref_start2 = pairs[i].alignment2.ref_start;
        int ref_id1 = pairs[i].alignment1.ref_id;
        int ref_id2 = pairs[i].alignment2.ref_id;
        if (
            ref_start1 != prev_ref_start1 ||
            ref_start2 != prev_ref_start2 ||
            ref_id1 != prev_ref_id1 ||
            ref_id2 != prev_ref_id2
        ) {
            prev_ref_start1 = ref_start1;
            prev_ref_start2 = ref_start2;
            prev_ref_id1 = ref_id1;
            prev_ref_id2 = ref_id2;
            pairs[j] = pairs[i];
            j++;
        }
    }
    pairs.resize(j);
}


/*
 * Count how many best alignments there are that all have the same score
 */
size_t count_best_alignment_pairs(const std::vector<ScoredAlignmentPair>& pairs) {
    if (pairs.empty()) {
        return 0;
    }
    size_t i = 1;
    for ( ; i < pairs.size(); ++i) {
        if (pairs[i].score != pairs[0].score) {
            break;
        }
    }
    return i;
}

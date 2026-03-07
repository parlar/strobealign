#ifndef STROBEALIGN_SV_TAGS_HPP
#define STROBEALIGN_SV_TAGS_HPP

#include <string>
#include <utility>
#include <vector>

#include "nam.hpp"
#include "refs.hpp"
#include "sam.hpp"
#include "statistics.hpp"

struct BreakpointInfo {
    int left_ref_id{-1}, right_ref_id{-1};
    int left_bp{0}, right_bp{0};
    int mh_fwd{0}, mh_bwd{0};
    int left_ci{0}, right_ci{0};
    std::string inserted_seq;
    std::string sv_type;

    std::string format_all(const References& refs) const;
};

std::string compute_yt_tag(
    const Alignment& a1, const Alignment& a2, float mu, float sigma);

std::string compute_yd_tag(
    const Alignment& a1, const Alignment& a2, float mu, float sigma);

std::string compute_xr_tag(const Details& details, bool is_fast_path);

std::string classify_sv_type(const Alignment& supp, const Alignment& primary);

std::string compute_ys_tag(const Alignment& supp, const Alignment& primary);

int facing_ref_pos(const Alignment& aln, bool covers_left_of_read);

std::pair<int, int> compute_microhomology(
    const References& references,
    int left_ref_id, int left_bp,
    int right_ref_id, int right_bp);

std::string compute_inserted_sequence(
    const Alignment& primary,
    const Alignment& supp,
    const std::string& read_seq);

int compute_ci_width(
    const References& references,
    int ref_id, int bp_pos,
    int microhomology_total);

BreakpointInfo compute_breakpoint_info(
    const Alignment& primary,
    const Alignment& supp,
    const References& references,
    const std::string& read_seq);

float compute_split_entropy(const std::vector<Alignment>& alignments);

float compute_max_artifact_prob(
    const Alignment& primary,
    const std::vector<Alignment>& supplementaries,
    int read_length);

std::string compute_supp_delta_tags(
    const Alignment& primary,
    const std::vector<Alignment>& supplementaries);

std::string build_sv_tags(
    const Alignment& a1, const Alignment& a2,
    float mu, float sigma,
    const Details& details, bool is_fast_path);

struct XuPosteriors {
    float ref;
    float sv;
    float rpt;
};

XuPosteriors compute_xu_posteriors(
    float score_ratio,
    float artifact_prob,
    float entropy,
    float mappability,
    int tied_count);

std::string format_xu_tag(const XuPosteriors& xu);

/*
 * Classify the mapping region repeat type from NAM distribution and hit stats.
 * Returns: 'U' (unique), 'T' (tandem), 'S' (segdup), 'D' (dispersed).
 */
char classify_repeat_type(
    const std::vector<Nam>& nams,
    const HitsDetails& hits);

#endif

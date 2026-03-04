#ifndef STROBEALIGN_CLIP_HPP
#define STROBEALIGN_CLIP_HPP

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "aligner.hpp"
#include "refs.hpp"
#include "sam.hpp"

constexpr int MIN_CLIP_LEN = 5;
constexpr int MAX_CLIP_LEN = 30;
constexpr int CLIP_SEARCH_WINDOW = 500;
constexpr float MIN_CLIP_IDENTITY = 0.75f;
constexpr int REFINE_PADDING = 10;

struct ClipRealignment {
    char side;        // '5' (leading) or '3' (trailing)
    int ref_pos;      // 1-based reference position
    std::string cigar;
    float identity;   // matches / alignment_length
};

/*
 * Slide clip across ref_window, count mismatches at each position.
 * Returns (best_position, mismatch_count) or (-1, -1) if clip is longer than window.
 */
std::pair<int, int> hamming_scan(const std::string& clip, const std::string& ref_window);

/*
 * Try realigning one clip against a reference window using Hamming scan + SSW refinement.
 * Returns a ClipRealignment if identity >= MIN_CLIP_IDENTITY.
 */
std::optional<ClipRealignment> try_realign_clip(
    const std::string& clip_seq,
    const std::string& ref_window,
    int window_offset,  // 0-based chromosome offset of window start
    char side,
    const Aligner& aligner);

/*
 * Realign short soft clips (5–30bp) against nearby reference.
 * Skips clip ends already covered by supplementary alignments.
 * read_seq must be on the mapped strand.
 */
std::vector<ClipRealignment> realign_soft_clips(
    const Alignment& primary,
    const std::string& read_seq,
    const References& references,
    const Aligner& aligner,
    bool left_covered = false,
    bool right_covered = false);

/* Format clip realignments as XK:Z tag value. */
std::string format_xk_tag(const std::vector<ClipRealignment>& clips);

#endif

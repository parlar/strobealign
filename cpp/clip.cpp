#include "clip.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cmath>

std::pair<int, int> hamming_scan(const std::string& clip, const std::string& ref_window) {
    int clen = static_cast<int>(clip.size());
    int wlen = static_cast<int>(ref_window.size());
    if (clen == 0 || clen > wlen) {
        return {-1, -1};
    }

    int best_pos = 0;
    int best_mm = clen;
    for (int pos = 0; pos <= wlen - clen; pos++) {
        int mm = 0;
        for (int i = 0; i < clen; i++) {
            if (std::toupper(clip[i]) != std::toupper(ref_window[pos + i])) {
                mm++;
            }
        }
        if (mm < best_mm) {
            best_mm = mm;
            best_pos = pos;
            if (mm == 0) break;
        }
    }
    return {best_pos, best_mm};
}

std::optional<ClipRealignment> try_realign_clip(
    const std::string& clip_seq,
    const std::string& ref_window,
    int window_offset,
    char side,
    const Aligner& aligner
) {
    auto [best_pos, best_mm] = hamming_scan(clip_seq, ref_window);
    if (best_pos < 0) return std::nullopt;

    int clip_len = static_cast<int>(clip_seq.size());
    float hamming_identity = static_cast<float>(clip_len - best_mm) / clip_len;

    // Too divergent even with indels — skip SSW
    if (hamming_identity < MIN_CLIP_IDENTITY * 0.8f) {
        return std::nullopt;
    }

    // Perfect or near-perfect Hamming match: build simple CIGAR
    if (hamming_identity >= MIN_CLIP_IDENTITY && best_mm <= 1) {
        // Build M CIGAR (simple match/mismatch)
        std::string cigar;
        if (best_mm == 0) {
            cigar = std::to_string(clip_len) + "M";
        } else {
            // Build EQ/X CIGAR showing exact match positions
            // For simplicity, just use nM
            cigar = std::to_string(clip_len) + "M";
        }
        return ClipRealignment{
            side,
            window_offset + best_pos + 1,  // 1-based
            cigar,
            hamming_identity
        };
    }

    // SSW refinement on narrow window around Hamming-best position
    int refine_start = std::max(0, best_pos - REFINE_PADDING);
    int refine_end = std::min(static_cast<int>(ref_window.size()), best_pos + clip_len + REFINE_PADDING);
    if (refine_end - refine_start < clip_len) {
        return std::nullopt;
    }

    std::string refine_window = ref_window.substr(refine_start, refine_end - refine_start);
    auto aln_result = aligner.align(clip_seq, refine_window);
    if (!aln_result) return std::nullopt;

    // Compute identity from alignment
    int matches = 0;
    int aln_len = 0;
    // Initialize q_pos=0 (not query_start) because the CIGAR already includes
    // a leading soft clip of length query_start when query_start > 0.
    int q_pos = 0;
    int r_pos = aln_result->ref_start;
    for (auto packed : aln_result->cigar.m_ops) {
        auto op = packed & 0xf;
        auto len = static_cast<int>(packed >> 4);
        switch (op) {
            case CIGAR_EQ:
                matches += len;
                aln_len += len;
                q_pos += len;
                r_pos += len;
                break;
            case CIGAR_X:
                aln_len += len;
                q_pos += len;
                r_pos += len;
                break;
            case CIGAR_MATCH:
                // M op: count actual matches by comparing sequences
                for (int i = 0; i < len; i++) {
                    if (q_pos + i < static_cast<int>(clip_seq.size())
                        && r_pos + i < static_cast<int>(refine_window.size())
                        && std::toupper(clip_seq[q_pos + i]) == std::toupper(refine_window[r_pos + i])) {
                        matches++;
                    }
                }
                aln_len += len;
                q_pos += len;
                r_pos += len;
                break;
            case CIGAR_INS:
                aln_len += len;
                q_pos += len;
                break;
            case CIGAR_DEL:
                aln_len += len;
                r_pos += len;
                break;
            case CIGAR_SOFTCLIP:
                q_pos += len;
                break;
            default:
                break;
        }
    }

    if (aln_len == 0) return std::nullopt;
    float identity = static_cast<float>(matches) / aln_len;
    if (identity < MIN_CLIP_IDENTITY) return std::nullopt;

    int mapped_pos = window_offset + refine_start + static_cast<int>(aln_result->ref_start) + 1;
    return ClipRealignment{
        side,
        mapped_pos,
        aln_result->cigar.to_string(),
        identity
    };
}

std::vector<ClipRealignment> realign_soft_clips(
    const Alignment& primary,
    const std::string& read_seq,
    const References& references,
    const Aligner& aligner,
    bool left_covered,
    bool right_covered
) {
    std::vector<ClipRealignment> results;
    if (primary.is_unaligned || primary.ref_id < 0) return results;

    int read_len = static_cast<int>(read_seq.size());
    // query_start/query_end are forward-strand coordinates, but read_seq is on
    // the mapped strand. For revcomp reads, forward-strand left = mapped-strand
    // right and vice versa, so swap the clip lengths.
    int left_clip, right_clip;
    if (primary.is_revcomp) {
        left_clip = read_len - primary.query_end;
        right_clip = primary.query_start;
    } else {
        left_clip = primary.query_start;
        right_clip = read_len - primary.query_end;
    }
    const auto& ref_seq = references.sequences[primary.ref_id];
    int ref_len = static_cast<int>(ref_seq.size());

    // Leading clip (5' in reference coordinates)
    if (!left_covered && left_clip >= MIN_CLIP_LEN && left_clip <= MAX_CLIP_LEN) {
        std::string clip_seq = read_seq.substr(0, left_clip);
        int window_end = primary.ref_start;
        int window_start = std::max(0, window_end - CLIP_SEARCH_WINDOW);
        if (window_end > window_start) {
            std::string ref_window = ref_seq.substr(window_start, window_end - window_start);
            auto cr = try_realign_clip(clip_seq, ref_window, window_start, '5', aligner);
            if (cr) results.push_back(*cr);
        }
    }

    // Trailing clip (3' in reference coordinates)
    if (!right_covered && right_clip >= MIN_CLIP_LEN && right_clip <= MAX_CLIP_LEN) {
        std::string clip_seq = read_seq.substr(read_len - right_clip);
        int aln_end = primary.ref_start + primary.length;
        int window_end = std::min(ref_len, aln_end + CLIP_SEARCH_WINDOW);
        if (window_end > aln_end) {
            std::string ref_window = ref_seq.substr(aln_end, window_end - aln_end);
            auto cr = try_realign_clip(clip_seq, ref_window, aln_end, '3', aligner);
            if (cr) results.push_back(*cr);
        }
    }

    return results;
}

std::string format_xk_tag(const std::vector<ClipRealignment>& clips) {
    std::string result;
    for (size_t i = 0; i < clips.size(); i++) {
        if (i > 0) result += ';';
        char buf[16];
        std::snprintf(buf, sizeof(buf), "%.2f", clips[i].identity);
        result += clips[i].side;
        result += ',';
        result += std::to_string(clips[i].ref_pos);
        result += ',';
        result += clips[i].cigar;
        result += ',';
        result += buf;
    }
    return result;
}

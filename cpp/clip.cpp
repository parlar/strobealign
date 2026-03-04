#include "clip.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cmath>
#include <cstring>

#ifdef __x86_64__
#include <immintrin.h>
#include <cpuid.h>
#endif

// Case-insensitive byte: uppercase ASCII letters via bitwise AND
static inline char upper_byte(char c) {
    return c & 0xDF;
}

static void hamming_scan_scalar(
    const char* clip_upper, const char* ref_upper,
    int clen, int wlen, int& best_pos, int& best_mm
) {
    for (int pos = 0; pos <= wlen - clen; pos++) {
        int mm = 0;
        for (int i = 0; i < clen; i++) {
            if (clip_upper[i] != ref_upper[pos + i]) {
                mm++;
            }
        }
        if (mm < best_mm) {
            best_mm = mm;
            best_pos = pos;
            if (mm == 0) break;
        }
    }
}

#ifdef __x86_64__
__attribute__((target("avx2")))
static void hamming_scan_avx2(
    const char* clip_upper, const char* ref_upper,
    int clen, int wlen, int& best_pos, int& best_mm
) {
    alignas(32) char clip_padded[32] = {};
    std::memcpy(clip_padded, clip_upper, clen);
    __m256i v_clip = _mm256_load_si256(reinterpret_cast<const __m256i*>(clip_padded));

    alignas(32) char mask_bytes[32] = {};
    std::memset(mask_bytes, 0xFF, clen);
    __m256i v_mask = _mm256_load_si256(reinterpret_cast<const __m256i*>(mask_bytes));

    for (int pos = 0; pos <= wlen - clen; pos++) {
        __m256i v_ref = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(ref_upper + pos));
        __m256i v_cmp = _mm256_cmpeq_epi8(v_clip, v_ref);
        v_cmp = _mm256_or_si256(v_cmp, _mm256_andnot_si256(v_mask, _mm256_set1_epi8(-1)));
        unsigned match_mask = static_cast<unsigned>(_mm256_movemask_epi8(v_cmp));
        int mm = 32 - __builtin_popcount(match_mask);

        if (mm < best_mm) {
            best_mm = mm;
            best_pos = pos;
            if (mm == 0) break;
        }
    }
}

__attribute__((target("sse2")))
static void hamming_scan_sse2(
    const char* clip_upper, const char* ref_upper,
    int clen, int wlen, int& best_pos, int& best_mm
) {
    alignas(16) char clip_padded[16] = {};
    std::memcpy(clip_padded, clip_upper, clen);
    __m128i v_clip = _mm_load_si128(reinterpret_cast<const __m128i*>(clip_padded));

    alignas(16) char mask_bytes[16] = {};
    std::memset(mask_bytes, 0xFF, clen);
    __m128i v_mask = _mm_load_si128(reinterpret_cast<const __m128i*>(mask_bytes));

    for (int pos = 0; pos <= wlen - clen; pos++) {
        __m128i v_ref = _mm_loadu_si128(reinterpret_cast<const __m128i*>(ref_upper + pos));
        __m128i v_cmp = _mm_cmpeq_epi8(v_clip, v_ref);
        v_cmp = _mm_or_si128(v_cmp, _mm_andnot_si128(v_mask, _mm_set1_epi8(-1)));
        int match_mask = _mm_movemask_epi8(v_cmp);
        int matches = __builtin_popcount(match_mask & 0xFFFF);
        int mm = 16 - matches;

        if (mm < best_mm) {
            best_mm = mm;
            best_pos = pos;
            if (mm == 0) break;
        }
    }
}

static bool cpu_has_avx2() {
    unsigned eax, ebx, ecx, edx;
    if (__get_cpuid_count(7, 0, &eax, &ebx, &ecx, &edx)) {
        return (ebx >> 5) & 1;  // AVX2 bit
    }
    return false;
}
#endif

std::pair<int, int> hamming_scan(const std::string& clip, const std::string& ref_window) {
    int clen = static_cast<int>(clip.size());
    int wlen = static_cast<int>(ref_window.size());
    if (clen == 0 || clen > wlen) {
        return {-1, -1};
    }

    // Pre-uppercase the clip (small, max 30bp)
    alignas(32) char clip_upper[64];
    for (int i = 0; i < clen; i++) {
        clip_upper[i] = upper_byte(clip[i]);
    }

    // Pre-uppercase the ref window (pad to allow 32-byte SIMD reads past end)
    std::string ref_upper(wlen + 32, '\0');
    for (int i = 0; i < wlen; i++) {
        ref_upper[i] = upper_byte(ref_window[i]);
    }

    int best_pos = 0;
    int best_mm = clen;

#ifdef __x86_64__
    static const bool has_avx2 = cpu_has_avx2();
    if (has_avx2) {
        hamming_scan_avx2(clip_upper, ref_upper.data(), clen, wlen, best_pos, best_mm);
    } else if (clen <= 16) {
        hamming_scan_sse2(clip_upper, ref_upper.data(), clen, wlen, best_pos, best_mm);
    } else {
        hamming_scan_scalar(clip_upper, ref_upper.data(), clen, wlen, best_pos, best_mm);
    }
#else
    hamming_scan_scalar(clip_upper, ref_upper.data(), clen, wlen, best_pos, best_mm);
#endif

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
        // clip_seq is a prefix of read_seq, must be a string for SSW compatibility
        std::string clip_seq(read_seq, 0, left_clip);
        int window_end = primary.ref_start;
        int window_start = std::max(0, window_end - CLIP_SEARCH_WINDOW);
        if (window_end > window_start) {
            // ref_window needs to be contiguous for SSW, use substr
            std::string ref_window(ref_seq, window_start, window_end - window_start);
            auto cr = try_realign_clip(clip_seq, ref_window, window_start, '5', aligner);
            if (cr) results.push_back(*cr);
        }
    }

    // Trailing clip (3' in reference coordinates)
    if (!right_covered && right_clip >= MIN_CLIP_LEN && right_clip <= MAX_CLIP_LEN) {
        std::string clip_seq(read_seq, read_len - right_clip);
        int aln_end = primary.ref_start + primary.length;
        int window_end = std::min(ref_len, aln_end + CLIP_SEARCH_WINDOW);
        if (window_end > aln_end) {
            std::string ref_window(ref_seq, aln_end, window_end - aln_end);
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

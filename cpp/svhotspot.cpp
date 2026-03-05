#include "svhotspot.hpp"
#include <algorithm>
#include <cmath>

void HotspotMap::build(
    const std::vector<SvEvidence>& evidence,
    int n_refs,
    const std::vector<unsigned int>& ref_lengths,
    int window_size,
    int min_support
) {
    intervals_.clear();
    intervals_.resize(n_refs);

    if (evidence.empty()) return;

    // Sort by (ref_id, position)
    auto sorted = evidence;
    std::sort(sorted.begin(), sorted.end(), [](const SvEvidence& a, const SvEvidence& b) {
        return a.ref_id < b.ref_id || (a.ref_id == b.ref_id && a.position < b.position);
    });

    // Process each reference separately
    size_t i = 0;
    while (i < sorted.size()) {
        int ref_id = sorted[i].ref_id;
        if (ref_id < 0 || ref_id >= n_refs) {
            i++;
            continue;
        }

        // Find range for this ref_id
        size_t j = i;
        while (j < sorted.size() && sorted[j].ref_id == ref_id) {
            j++;
        }

        // Scale min_support based on evidence density for this reference.
        // If background evidence rate is high, raise the threshold to only
        // flag windows significantly above the background.
        size_t ref_evidence_count = j - i;
        int effective_min_support = min_support;
        if (ref_id < static_cast<int>(ref_lengths.size()) && ref_lengths[ref_id] > 0) {
            double density = static_cast<double>(ref_evidence_count) / ref_lengths[ref_id];
            double expected_per_window = density * window_size;
            // Require at least 5x the expected background rate
            int density_threshold = static_cast<int>(std::ceil(expected_per_window * 5.0));
            effective_min_support = std::max(min_support, density_threshold);
        }

        // Sliding window: for each evidence item, count how many items
        // fall within [position, position + window_size)
        auto& ref_intervals = intervals_[ref_id];
        size_t left = i;
        for (size_t right = i; right < j; right++) {
            // Advance left pointer past the window
            while (sorted[left].position < sorted[right].position - window_size) {
                left++;
            }
            int count = static_cast<int>(right - left + 1);
            if (count >= effective_min_support) {
                int start = sorted[left].position;
                int end = sorted[right].position;

                // Merge with previous interval if overlapping/adjacent
                if (!ref_intervals.empty() && start <= ref_intervals.back().end + window_size) {
                    ref_intervals.back().end = std::max(ref_intervals.back().end, end);
                    ref_intervals.back().support = std::max(ref_intervals.back().support, count);
                } else {
                    ref_intervals.push_back({start, end, count});
                }
            }
        }

        i = j;
    }
}

bool HotspotMap::near_hotspot(int ref_id, int position, int margin) const {
    if (ref_id < 0 || ref_id >= static_cast<int>(intervals_.size())) {
        return false;
    }
    const auto& ref_intervals = intervals_[ref_id];
    if (ref_intervals.empty()) return false;

    // Binary search: find first interval whose end >= position - margin
    int query_start = position - margin;
    int query_end = position + margin;

    auto it = std::lower_bound(
        ref_intervals.begin(), ref_intervals.end(), query_start,
        [](const HotspotInterval& interval, int val) {
            return interval.end < val;
        }
    );

    // Check if this interval overlaps [query_start, query_end]
    if (it != ref_intervals.end() && it->start <= query_end) {
        return true;
    }
    return false;
}

size_t HotspotMap::total_hotspots() const {
    size_t total = 0;
    for (const auto& ref : intervals_) {
        total += ref.size();
    }
    return total;
}

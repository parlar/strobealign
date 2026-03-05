#ifndef STROBEALIGN_SVHOTSPOT_HPP
#define STROBEALIGN_SVHOTSPOT_HPP

#include <vector>
#include <cstddef>

struct SvEvidence {
    int ref_id;
    int position;
};

/*
 * Thread-local SV evidence collector. No synchronization needed —
 * each worker thread owns its own instance.
 */
struct SvEvidenceCollector {
    std::vector<SvEvidence> evidence;

    void add(int ref_id, int pos) {
        evidence.push_back({ref_id, pos});
    }
};

struct HotspotInterval {
    int start;
    int end;
    int support;
};

/*
 * Sorted collection of genomic intervals with high SV evidence density.
 * Built once between pass 1 and pass 2, then queried read-only during pass 2.
 */
class HotspotMap {
public:
    /*
     * Build hotspot intervals from merged evidence.
     * 1. Sort evidence by (ref_id, position).
     * 2. Slide a window of `window_size` bp, count evidence in each window.
     * 3. Windows with >= min_support items become hotspot intervals.
     * 4. Merge overlapping/adjacent intervals.
     */
    void build(const std::vector<SvEvidence>& evidence,
               int n_refs, const std::vector<unsigned int>& ref_lengths,
               int window_size = 500, int min_support = 3);

    /* Check if `position` on `ref_id` is within `margin` bp of any hotspot. */
    bool near_hotspot(int ref_id, int position, int margin = 500) const;

    size_t total_hotspots() const;

private:
    // Per-reference sorted intervals. Indexed by ref_id.
    std::vector<std::vector<HotspotInterval>> intervals_;
};

#endif

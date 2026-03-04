#ifndef STROBEALIGN_PHASING_HPP
#define STROBEALIGN_PHASING_HPP

#include <string>
#include <vector>
#include <optional>

#include "refs.hpp"
#include "sam.hpp"

struct PhasedVariant {
    int position;    // 0-based reference position
    char ref_allele; // Reference base (uppercase)
    char alt_allele; // Alternate base (uppercase)
    int phase;       // 1 if haplotype 1 carries alt, 2 if haplotype 2 carries alt
    int phase_set;   // Phase set ID (from VCF PS field or positional block)
};

struct PhaseResult {
    int haplotype;   // 1 or 2
    int phase_set;   // Phase set ID
    int support;     // Number of informative sites supporting this assignment
};

/*
 * Read-only map of phased heterozygous variants, indexed by reference ID.
 * Built once from a phased VCF, then queried concurrently by worker threads.
 */
class PhasingMap {
public:
    /*
     * Load phased heterozygous SNPs from a VCF file.
     * Only loads biallelic SNPs with phased GT (pipe-separated: 0|1, 1|0).
     * Requires refs for chromosome name -> ref_id mapping.
     */
    void load(const std::string& vcf_filename, const References& refs);

    /*
     * Assign haplotype to a read based on overlapping phased variants.
     * Returns nullopt if no informative sites overlap or if evidence is tied.
     */
    std::optional<PhaseResult> assign_haplotype(
        const Alignment& aln,
        const std::string& read_seq,
        const Cigar& cigar) const;

    size_t total_variants() const;

    bool empty() const { return variants_.empty(); }

private:
    // Per-reference sorted variants. Indexed by ref_id.
    std::vector<std::vector<PhasedVariant>> variants_;
};

std::string format_phase_tags(const PhaseResult& result);

#endif

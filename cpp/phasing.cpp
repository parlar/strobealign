#include "phasing.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string_view>
#include <unordered_map>

#include "logger.hpp"

static Logger& logger = Logger::get();

// Fast tab-field extraction: returns string_view to the n-th tab-separated field (0-based)
// Returns empty view if field not found
static inline std::string_view get_tab_field(const char* line, size_t len, int field_idx) {
    const char* p = line;
    const char* end = line + len;
    for (int i = 0; i < field_idx; i++) {
        p = static_cast<const char*>(std::memchr(p, '\t', end - p));
        if (!p) return {};
        p++;
    }
    const char* tab = static_cast<const char*>(std::memchr(p, '\t', end - p));
    return {p, static_cast<size_t>(tab ? tab - p : end - p)};
}

// Fast colon-field extraction: returns string_view to the n-th colon-separated field
static inline std::string_view get_colon_field(std::string_view sv, int field_idx) {
    const char* p = sv.data();
    const char* end = p + sv.size();
    for (int i = 0; i < field_idx; i++) {
        p = static_cast<const char*>(std::memchr(p, ':', end - p));
        if (!p) return {};
        p++;
    }
    const char* colon = static_cast<const char*>(std::memchr(p, ':', end - p));
    return {p, static_cast<size_t>(colon ? colon - p : end - p)};
}

// Find the index of a field name (like "GT" or "PS") in a colon-separated FORMAT string
static inline int find_format_index(std::string_view format, std::string_view name) {
    const char* p = format.data();
    const char* end = p + format.size();
    int idx = 0;
    while (p < end) {
        const char* colon = static_cast<const char*>(std::memchr(p, ':', end - p));
        size_t flen = colon ? static_cast<size_t>(colon - p) : static_cast<size_t>(end - p);
        if (flen == name.size() && std::memcmp(p, name.data(), flen) == 0) {
            return idx;
        }
        if (!colon) break;
        p = colon + 1;
        idx++;
    }
    return -1;
}

/*
 * Parse a phased VCF file and load biallelic SNPs with phased genotypes.
 * Expected GT format: 0|1 or 1|0 (pipe = phased).
 * Skips: unphased (0/1), homozygous (0|0, 1|1), multiallelic, indels.
 */
void PhasingMap::load(const std::string& vcf_filename, const References& refs) {
    // Build name -> ref_id lookup
    std::unordered_map<std::string, int> name_to_id;
    for (size_t i = 0; i < refs.names.size(); ++i) {
        name_to_id[refs.names[i]] = static_cast<int>(i);
    }

    variants_.resize(refs.names.size());

    std::ifstream vcf(vcf_filename);
    if (!vcf.is_open()) {
        throw std::runtime_error("Cannot open VCF file: " + vcf_filename);
    }

    size_t loaded = 0;
    size_t skipped = 0;
    std::string line;
    while (std::getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;

        const char* data = line.data();
        size_t len = line.size();

        // Fields: CHROM(0) POS(1) ID(2) REF(3) ALT(4) QUAL(5) FILTER(6) INFO(7) FORMAT(8) SAMPLE(9)
        auto chrom_sv = get_tab_field(data, len, 0);
        auto pos_sv = get_tab_field(data, len, 1);
        auto ref_sv = get_tab_field(data, len, 3);
        auto alt_sv = get_tab_field(data, len, 4);
        auto format_sv = get_tab_field(data, len, 8);
        auto sample_sv = get_tab_field(data, len, 9);

        if (chrom_sv.empty() || pos_sv.empty() || ref_sv.empty() || alt_sv.empty()
            || format_sv.empty() || sample_sv.empty()) {
            skipped++;
            continue;
        }

        // Skip non-SNPs
        if (ref_sv.size() != 1 || alt_sv.size() != 1 || alt_sv.find(',') != std::string_view::npos) {
            skipped++;
            continue;
        }

        // Parse position
        int pos = 0;
        auto [ptr, ec] = std::from_chars(pos_sv.data(), pos_sv.data() + pos_sv.size(), pos);
        if (ec != std::errc()) {
            skipped++;
            continue;
        }

        // Find GT and PS indices in FORMAT
        int gt_index = find_format_index(format_sv, "GT");
        if (gt_index < 0) {
            skipped++;
            continue;
        }
        int ps_index = find_format_index(format_sv, "PS");

        // Extract GT value
        auto gt_val = get_colon_field(sample_sv, gt_index);

        // Must be phased (pipe-separated) heterozygous
        if (gt_val.size() != 3 || gt_val[1] != '|') {
            skipped++;
            continue;
        }
        int allele1 = gt_val[0] - '0';
        int allele2 = gt_val[2] - '0';
        if (allele1 == allele2) {
            skipped++;
            continue;
        }
        if ((allele1 != 0 && allele1 != 1) || (allele2 != 0 && allele2 != 1)) {
            skipped++;
            continue;
        }

        // Extract PS value if available
        int ps_val = 0;
        if (ps_index >= 0) {
            auto ps_sv = get_colon_field(sample_sv, ps_index);
            if (!ps_sv.empty()) {
                std::from_chars(ps_sv.data(), ps_sv.data() + ps_sv.size(), ps_val);
            }
        }

        // Find ref_id
        std::string chrom_str(chrom_sv);
        auto it = name_to_id.find(chrom_str);
        if (it == name_to_id.end()) {
            skipped++;
            continue;
        }

        int phase = (allele1 == 1) ? 1 : 2;
        if (ps_val == 0) ps_val = pos;

        PhasedVariant var;
        var.position = pos - 1;  // VCF is 1-based, convert to 0-based
        var.ref_allele = ref_sv[0] & 0xDF;  // uppercase via bitmask
        var.alt_allele = alt_sv[0] & 0xDF;
        var.phase = phase;
        var.phase_set = ps_val;

        variants_[it->second].push_back(var);
        loaded++;
    }

    // Sort each chromosome's variants by position
    for (auto& vars : variants_) {
        std::sort(vars.begin(), vars.end(),
            [](const PhasedVariant& a, const PhasedVariant& b) {
                return a.position < b.position;
            });
    }

    logger.debug() << "Phasing: loaded " << loaded << " het SNPs, skipped " << skipped << "\n";
}

/*
 * Walk the CIGAR to find which base the read has at each phased variant position.
 * Count votes for haplotype 1 vs haplotype 2.
 */
static inline const char* make_complement_table() {
    static char table[256] = {};
    table[static_cast<unsigned char>('A')] = 'T';
    table[static_cast<unsigned char>('T')] = 'A';
    table[static_cast<unsigned char>('C')] = 'G';
    table[static_cast<unsigned char>('G')] = 'C';
    table[static_cast<unsigned char>('a')] = 't';
    table[static_cast<unsigned char>('t')] = 'a';
    table[static_cast<unsigned char>('c')] = 'g';
    table[static_cast<unsigned char>('g')] = 'c';
    return table;
}

static inline char complement(char c) {
    static const char* table = make_complement_table();
    char r = table[static_cast<unsigned char>(c)];
    return r ? r : 'N';
}

static std::string reverse_complement(const std::string& seq) {
    std::string rc(seq.size(), 'N');
    for (size_t i = 0; i < seq.size(); ++i) {
        rc[seq.size() - 1 - i] = complement(seq[i]);
    }
    return rc;
}

std::optional<PhaseResult> PhasingMap::assign_haplotype(
    const Alignment& aln,
    const std::string& read_seq,
    const Cigar& cigar
) const {
    if (aln.is_unaligned || aln.ref_id < 0) return std::nullopt;
    if (static_cast<size_t>(aln.ref_id) >= variants_.size()) return std::nullopt;
    const auto& vars = variants_[aln.ref_id];
    if (vars.empty()) return std::nullopt;

    // CIGAR is relative to the mapped strand. If the read mapped to the reverse
    // strand, we need to RC the original FASTQ sequence to get the aligned bases.
    const std::string rc_seq = aln.is_revcomp ? reverse_complement(read_seq) : std::string();
    const std::string& query = aln.is_revcomp ? rc_seq : read_seq;

    int ref_end = aln.ref_start + aln.length;

    // Binary search for first variant at or after ref_start
    auto it = std::lower_bound(vars.begin(), vars.end(), aln.ref_start,
        [](const PhasedVariant& v, int pos) { return v.position < pos; });

    if (it == vars.end() || it->position >= ref_end) return std::nullopt;

    // Walk CIGAR, tracking ref and query positions
    int r_pos = aln.ref_start;
    int q_pos = 0;
    int h1_votes = 0, h2_votes = 0;
    int phase_set = 0;

    for (auto packed : cigar.m_ops) {
        auto op = packed & 0xf;
        auto len = static_cast<int>(packed >> 4);

        switch (op) {
            case CIGAR_MATCH:
            case CIGAR_EQ:
            case CIGAR_X:
                // Check variants in this M/=/X block
                while (it != vars.end() && it->position < r_pos + len) {
                    if (it->position >= r_pos && it->position < ref_end) {
                        int offset = it->position - r_pos;
                        if (q_pos + offset >= 0 && q_pos + offset < static_cast<int>(query.size())) {
                            char read_base = std::toupper(query[q_pos + offset]);
                            if (read_base == it->alt_allele) {
                                // Read supports the ALT allele
                                if (it->phase == 1) h1_votes++;
                                else h2_votes++;
                            } else if (read_base == it->ref_allele) {
                                // Read supports the REF allele
                                if (it->phase == 1) h2_votes++;
                                else h1_votes++;
                            }
                            if (phase_set == 0) phase_set = it->phase_set;
                        }
                    }
                    ++it;
                }
                r_pos += len;
                q_pos += len;
                break;
            case CIGAR_INS:
                q_pos += len;
                break;
            case CIGAR_DEL:
                // Skip variants within deletion
                while (it != vars.end() && it->position < r_pos + len) {
                    ++it;
                }
                r_pos += len;
                break;
            case CIGAR_SOFTCLIP:
                q_pos += len;
                break;
            default:
                break;
        }
    }

    if (h1_votes == 0 && h2_votes == 0) return std::nullopt;
    if (h1_votes == h2_votes) return std::nullopt;  // Tied — can't assign

    int hap = (h1_votes > h2_votes) ? 1 : 2;
    int support = std::max(h1_votes, h2_votes);
    return PhaseResult{hap, phase_set, support};
}

size_t PhasingMap::total_variants() const {
    size_t total = 0;
    for (const auto& vars : variants_) {
        total += vars.size();
    }
    return total;
}

std::string format_phase_tags(const PhaseResult& result) {
    return "HP:i:" + std::to_string(result.haplotype)
         + "\tPS:i:" + std::to_string(result.phase_set);
}

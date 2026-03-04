#include "phasing.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#include "logger.hpp"

static Logger& logger = Logger::get();

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

        // Parse tab-separated fields: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE...
        std::istringstream iss(line);
        std::string chrom, id, ref, alt, qual, filter, info, format, sample;
        int pos;
        if (!(iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format >> sample)) {
            skipped++;
            continue;
        }

        // Skip non-SNPs (indels, multiallelic)
        if (ref.size() != 1 || alt.size() != 1 || alt.find(',') != std::string::npos) {
            skipped++;
            continue;
        }

        // Find GT field index in FORMAT
        int gt_index = -1;
        int ps_index = -1;
        {
            std::istringstream fmt_iss(format);
            std::string field;
            int idx = 0;
            while (std::getline(fmt_iss, field, ':')) {
                if (field == "GT") gt_index = idx;
                if (field == "PS") ps_index = idx;
                idx++;
            }
        }
        if (gt_index < 0) {
            skipped++;
            continue;
        }

        // Extract GT and PS from sample field
        std::string gt_val;
        int ps_val = 0;
        {
            std::istringstream smp_iss(sample);
            std::string field;
            int idx = 0;
            while (std::getline(smp_iss, field, ':')) {
                if (idx == gt_index) gt_val = field;
                if (idx == ps_index) {
                    try { ps_val = std::stoi(field); } catch (...) { ps_val = 0; }
                }
                idx++;
            }
        }

        // Must be phased (pipe-separated) heterozygous
        if (gt_val.size() != 3 || gt_val[1] != '|') {
            skipped++;
            continue;
        }
        int allele1 = gt_val[0] - '0';
        int allele2 = gt_val[2] - '0';
        if (allele1 == allele2) {
            skipped++;  // Homozygous
            continue;
        }
        if ((allele1 != 0 && allele1 != 1) || (allele2 != 0 && allele2 != 1)) {
            skipped++;  // Non-biallelic GT
            continue;
        }

        // Find ref_id
        auto it = name_to_id.find(chrom);
        if (it == name_to_id.end()) {
            skipped++;
            continue;
        }

        // phase: which haplotype carries the ALT allele
        // If GT=0|1, haplotype 2 carries ALT
        // If GT=1|0, haplotype 1 carries ALT
        int phase = (allele1 == 1) ? 1 : 2;

        // Use PS field if available, otherwise use position as phase set
        if (ps_val == 0) ps_val = pos;

        PhasedVariant var;
        var.position = pos - 1;  // VCF is 1-based, convert to 0-based
        var.ref_allele = std::toupper(ref[0]);
        var.alt_allele = std::toupper(alt[0]);
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
static char complement(char c) {
    switch (c) {
        case 'A': return 'T'; case 'T': return 'A';
        case 'C': return 'G'; case 'G': return 'C';
        case 'a': return 't'; case 't': return 'a';
        case 'c': return 'g'; case 'g': return 'c';
        default: return 'N';
    }
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

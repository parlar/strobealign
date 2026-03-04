#ifndef STROBEALIGN_MAPPINGPARAMETERS_HPP
#define STROBEALIGN_MAPPINGPARAMETERS_HPP

#include "mcsstrategy.hpp"
#include "sam.hpp"
#include "exceptions.hpp"

enum class OutputFormat {
    SAM,
    PAF,
    Abundance
};

struct ChainingParameters {
    int max_lookback;
    float diag_diff_penalty;
    float gap_length_penalty;
    float valid_score_threshold;
    int max_ref_gap;
    float matches_weight;
};

struct MappingParameters {
    int r { 150 };
    int max_secondary { 0 };
    int max_supplementary { 0 };
    int max_supp_overlap { 50 };
    int min_clip { 15 };
    float dropoff_threshold { 0.5 };
    int rescue_threshold{100};
    int max_tries { 20 };
    McsStrategy mcs_strategy{McsStrategy::Always};
    OutputFormat output_format {OutputFormat::SAM};
    CigarOps cigar_ops{CigarOps::M};
    bool output_unmapped { true };
    bool details{false};
    bool fastq_comments{false};

    bool sv_tags{false};
    bool use_nams{false};
    ChainingParameters chaining_params;

    void verify() const {
        if (max_tries < 1) {
            throw BadParameter("max_tries must be greater than zero");
        }
        if (max_supplementary < 0) {
            throw BadParameter("--supp must be non-negative");
        }
        if (max_supp_overlap < 0) {
            throw BadParameter("--supp-overlap must be non-negative");
        }
        if (min_clip < 1) {
            throw BadParameter("--min-clip must be at least 1");
        }
    }
};

#endif

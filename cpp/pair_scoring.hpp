#ifndef STROBEALIGN_PAIR_SCORING_HPP
#define STROBEALIGN_PAIR_SCORING_HPP

#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

#include "aligner.hpp"
#include "nam.hpp"
#include "revcomp.hpp"
#include "sam.hpp"

struct NamPair {
    float score;
    Nam nam1;
    Nam nam2;
};

std::ostream& operator<<(std::ostream& os, const NamPair& nam_pair);

struct ScoredAlignmentPair {
    double score;
    Alignment alignment1;
    Alignment alignment2;
};

float compute_pair_posterior(const std::vector<ScoredAlignmentPair>& pairs);

uint8_t proper_pair_mapq(const std::vector<Nam>& nams);

std::pair<int, int> joint_mapq_from_high_scores(const std::vector<ScoredAlignmentPair>& pairs);

float normal_pdf(float x, float mu, float sigma);

std::vector<ScoredAlignmentPair> get_best_scoring_pairs(
    const std::vector<Alignment>& alignments1,
    const std::vector<Alignment>& alignments2,
    float mu,
    float sigma);

bool is_proper_nam_pair(const Nam nam1, const Nam nam2, float mu, float sigma);

std::vector<NamPair> get_best_scoring_nam_pairs(
    const std::vector<Nam>& nams1,
    const std::vector<Nam>& nams2,
    float mu,
    float sigma);

Alignment rescue_align(
    const Aligner& aligner,
    const Nam& mate_nam,
    const References& references,
    const Read& read,
    float mu,
    float sigma,
    int k);

void deduplicate_scored_pairs(std::vector<ScoredAlignmentPair>& pairs);

size_t count_best_alignment_pairs(const std::vector<ScoredAlignmentPair>& pairs);

#endif

#ifndef SLICEABLE_SET_H
#define SLICEABLE_SET_H

#include <cstdint>
#include <vector>

#include "common.hpp"

bool is_subset(const sliceable_set_t& subset,
               const std::vector<sliceable_set_t>& combos);

sliceable_set_t unique_sliceable_set_naive(
    const sliceable_set_t& ss, const std::vector<edge_trans_t>& transformations,
    int32_t n);

sliceable_set_t unique_sliceable_set(
    const sliceable_set_t& ss, const std::vector<edge_trans_t>& transformations,
    int32_t n);

std::vector<sliceable_set_t> combine_usr_mss(
    const std::vector<sliceable_set_t>& usr,
    const std::vector<sliceable_set_t>& mss,
    const std::vector<edge_trans_t>& transformations, int32_t n);

std::vector<sliceable_set_t> usr_to_mss(
    const std::vector<sliceable_set_t>& usr,
    const std::vector<edge_trans_t>& transformations, int32_t n);

#endif
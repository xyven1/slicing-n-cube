#ifndef LOW_WEIGHT_H
#define LOW_WEIGHT_H

#include <cstdint>
#include <unordered_set>
#include <vector>

#include "common.hpp"

template <int32_t N>
sliceable_set_t<N> low_weight_halfspace_to_sliceable_set(
    const std::vector<int32_t>& normal, int32_t distance,
    const std::vector<edge_t>& edges);

template <int32_t N>
std::vector<sliceable_set_t<N>> compute_low_weight_sliceable_sets(
    int32_t min_weight, int32_t max_weight, const std::vector<edge_t>& edges);

template <int32_t N>
int32_t combine_low_weight_sliceable_sets(
    const std::vector<sliceable_set_t<N>>& sets,
    const std::vector<edge_trans_t>& transformations);

#endif

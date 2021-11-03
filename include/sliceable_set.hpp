#ifndef SLICEABLE_SET_H
#define SLICEABLE_SET_H

#include <cstdint>
#include <filesystem>
#include <vector>

#include "common.hpp"

template <int32_t N>
bool is_subset(const sliceable_set_t<N>& subset,
               const std::vector<sliceable_set_t<N>>& combos);

template <int32_t N>
sliceable_set_t<N> unique_sliceable_set_naive(
    const sliceable_set_t<N>& ss,
    const std::vector<edge_trans_t>& transformations);

template <int32_t N>
sliceable_set_t<N> unique_sliceable_set(
    const sliceable_set_t<N>& ss,
    const std::vector<edge_trans_t>& transformations);

template <int32_t N>
std::vector<sliceable_set_t<N>> combine_usr_mss(
    const std::vector<sliceable_set_t<N>>& usr,
    const std::vector<sliceable_set_t<N>>& mss,
    const std::vector<edge_trans_t>& transformations);

template <int32_t N>
bool combine_usr_mss_final(const std::vector<sliceable_set_t<N>>& usr,
                           const std::vector<sliceable_set_t<N>>& mss);

template <int32_t N>
bool combine_usr_mss_final_naive(const std::vector<sliceable_set_t<N>>& usr,
                                 const std::vector<sliceable_set_t<N>>& mss);

template <int32_t N>
std::vector<sliceable_set_t<N>> usr_to_mss(
    const std::vector<sliceable_set_t<N>>& usr,
    const std::vector<edge_trans_t>& transformations);

template <int32_t N>
int32_t get_leading_zeros(const sliceable_set_t<N>& ss);

template <int32_t N>
int32_t get_leading_ones(const sliceable_set_t<N>& ss);

template <int32_t N>
void write_to_file(const std::vector<sliceable_set_t<N>>& sets,
                   const std::filesystem::path& path);

template <int32_t N>
std::vector<sliceable_set_t<N>> read_from_file(
    const std::filesystem::path& path);

#endif
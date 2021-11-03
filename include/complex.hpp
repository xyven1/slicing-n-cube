#ifndef COMPLEX_H
#define COMPLEX_H

#include <cstdint>
#include <vector>

#include "common.hpp"

template <int32_t N>
std::vector<vertex_t> adjacent_vertices_of_complex(const complex_t<N>& complex);

template <int32_t N>
complex_t<N> unique_complex(const complex_t<N>& complex,
                            const std::vector<vertex_trans_t>& transformations);

template <int32_t N>
sliceable_set_t<N> complex_to_sliceable_set(const complex_t<N>& complex,
                                            const std::vector<edge_t>& edges);

template <int32_t N>
std::vector<sliceable_set_t<N>> complexes_to_usr(
    const std::vector<complex_t<N>>& complexes,
    const std::vector<edge_trans_t>& edge_transformations,
    const std::vector<edge_t>& edges);

template <int32_t N>
std::vector<sliceable_set_t<N>> complexes_to_mss(
    const std::vector<complex_t<N>>& complexes,
    const std::vector<vertex_trans_t>& vertex_transformations,
    const std::vector<edge_t>& edges);

template <int32_t N>
std::vector<complex_t<N>> compute_cut_complexes(
    const std::vector<vertex_trans_t>& transformations);

#endif

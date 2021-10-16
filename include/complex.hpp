#ifndef COMPLEX_H
#define COMPLEX_H

#include <cstdint>
#include <vector>

#include "common.hpp"

std::vector<vertex_t> adjacent_vertices_of_complex(const complex_t& complex,
                                                   int32_t n);

complex_t unique_complex(const complex_t& complex,
                         const std::vector<vertex_trans_t>& transformations,
                         int32_t n);

sliceable_set_t complex_to_sliceable_set(const complex_t& complex,
                                         const std::vector<edge_t>& edges,
                                         int32_t n);

std::vector<complex_t> compute_cut_complexes(
    const std::vector<vertex_trans_t>& transformations, int32_t n);

std::vector<complex_t> expand_complex(
    const std::vector<vertex_trans_t>& transformations,
    const complex_t& complex, int32_t n);

#endif

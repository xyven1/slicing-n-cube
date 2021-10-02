#ifndef COMPLEX_H
#define COMPLEX_H

#include <cstdint>
#include <vector>

#include "common.hpp"

std::vector<vertex_t> adjacent_vertices_of_complex(const complex_t& complex,
                                                   int32_t n);

complex_t unique_complex(const complex_t& complex,
                         const std::vector<symmetry_t>& symmetries, int32_t n);

std::vector<edge_t> complex_to_edges(const complex_t& complex, int32_t n);

std::vector<complex_t> compute_cut_complexes(int32_t n);

#endif

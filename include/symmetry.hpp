#ifndef SYMMETRY_H
#define SYMMETRY_H

#include <cstdint>
#include <vector>

#include "common.hpp"

std::vector<symmetry_t> compute_symmetries(int32_t n);

vertex_t transform_vertex(const symmetry_t& sym, vertex_t v, int32_t n);

vertex_t transform_vertex_inv(const symmetry_t& sym, vertex_t v, int32_t n);

complex_t transform_complex(const complex_t& complex, const symmetry_t& sym,
                            int32_t n);

complex_t transform_complex_and_min(const complex_t& complex,
                                    const symmetry_t& sym, int32_t n,
                                    const complex_t& min_complex);

std::vector<transformation_t> compute_vertex_transformations(
    const std::vector<symmetry_t>& symmetries, int32_t n);

std::vector<inversion_t> compute_vertex_inversions(
    const std::vector<symmetry_t>& symmetries, int32_t n);

#endif

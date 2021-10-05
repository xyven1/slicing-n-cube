#ifndef SYMMETRY_H
#define SYMMETRY_H

#include <cstdint>
#include <vector>

#include "common.hpp"

std::vector<symmetry_t> compute_symmetries(int32_t n);

vertex_t transform_vertex(const symmetry_t& sym, vertex_t v, int32_t n);

vertex_t transform_vertex_inv(const symmetry_t& sym, vertex_t v, int32_t n);

std::vector<transformation_t> compute_vertex_transformations(
    const std::vector<symmetry_t>& symmetries, int32_t n);

std::vector<inversion_t> compute_vertex_inversions(
    const std::vector<symmetry_t>& symmetries, int32_t n);

#endif

#ifndef SYMMETRY_H
#define SYMMETRY_H

#include <cstdint>
#include <vector>

#include "common.hpp"

std::vector<symmetry_t> compute_symmetries(int32_t n);

vertex_t transform_vertex(const symmetry_t& sym, vertex_t v, int32_t n);

vertex_t transform_vertex_inv(const symmetry_t& sym, vertex_t v, int32_t n);

std::vector<vertex_trans_t> compute_vertex_transformations(int32_t n);

std::vector<edge_trans_t> compute_edge_transformations(
    const std::vector<edge_t>& edges,
    const std::vector<vertex_trans_t>& vertex_transformations, int32_t n);

#endif

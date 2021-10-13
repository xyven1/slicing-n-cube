#include "symmetry.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

#include "common.hpp"
#include "edge.hpp"

std::vector<symmetry_t> compute_symmetries(int32_t n) {
  std::vector<symmetry_t> symmetries;
  symmetries.reserve(num_symmetries(n));

  std::vector<int32_t> permutation(n);
  for (int32_t i = 0; i < n; ++i) {
    permutation[i] = i;
  }
  do {
    // if the i-th least significant bit is set, negate the i-th coordinate
    for (uint32_t signs = 0; signs < (1u << n); ++signs) {
      symmetry_t symmetry;
      symmetry.reserve(n);
      for (int32_t i = 0; i < n; ++i) {
        const uint32_t sign = (signs >> i) & 1;
        const uint32_t sign_and_position = (sign << 31) | permutation[i];
        symmetry.push_back(sign_and_position);
      }
      symmetries.push_back(symmetry);
    }
  } while (std::next_permutation(permutation.begin(), permutation.end()));
  return symmetries;
}

vertex_t transform_vertex(const symmetry_t& sym, vertex_t v, int32_t n) {
  vertex_t transformation = 0;
  for (int32_t i = 0; i < n; ++i) {
    const uint32_t sign = (sym[i] >> 31) & 1;
    const uint32_t position = sym[i] & 0x7FFFFFFF;
    const uint32_t val_i = (v >> i) & 1;
    const uint32_t val_i_sign = val_i ^ sign;
    const uint32_t val_i_sign_and_position = val_i_sign << position;
    transformation |= val_i_sign_and_position;
  }
  return transformation;
}

vertex_t transform_vertex_inv(const symmetry_t& sym, vertex_t v, int32_t n) {
  vertex_t transformation = 0;
  for (int32_t i = 0; i < n; ++i) {
    const uint32_t sign = (sym[i] >> 31) & 1;
    const uint32_t position = sym[i] & 0x7FFFFFFF;
    const uint32_t val_i = (v >> position) & 1;
    const uint32_t val_i_sign = val_i ^ sign;
    const uint32_t val_i_sign_and_position = val_i_sign << i;
    transformation |= val_i_sign_and_position;
  }
  return transformation;
}

std::vector<vertex_trans_t> compute_vertex_transformations(
    const std::vector<symmetry_t>& symmetries, int32_t n) {
  std::vector<vertex_trans_t> transformations(symmetries.size());
  for (std::size_t i = 0; i < symmetries.size(); ++i) {
    transformations[i].reserve(num_vertices(n));
    for (vertex_t v = 0; v < num_vertices(n); ++v) {
      transformations[i].push_back(transform_vertex(symmetries[i], v, n));
    }
  }
  return transformations;
}

std::vector<vertex_inv_t> compute_vertex_inversions(
    const std::vector<symmetry_t>& symmetries, int32_t n) {
  std::vector<vertex_inv_t> inversions(symmetries.size());
  for (std::size_t i = 0; i < symmetries.size(); ++i) {
    inversions[i].reserve(num_vertices(n));
    for (vertex_t v = 0; v < num_vertices(n); ++v) {
      inversions[i].push_back(transform_vertex_inv(symmetries[i], v, n));
    }
  }
  return inversions;
}

std::vector<edge_trans_t> compute_edge_transformations(
    const std::vector<edge_t>& edges,
    const std::vector<vertex_trans_t>& vertex_transformations, int32_t n) {
  std::vector<edge_trans_t> edge_transformations(vertex_transformations.size());
  for (std::size_t i = 0; i < edge_transformations.size(); ++i) {
    edge_transformations[i].reserve(num_edges(n));
    for (int32_t e = 0; e < num_edges(n); ++e) {
      const vertex_t u_trans = vertex_transformations[i][edges[e].first];
      const vertex_t v_trans = vertex_transformations[i][edges[e].second];
      const edge_t e_trans = (u_trans < v_trans) ? edge_t(u_trans, v_trans)
                                                 : edge_t(v_trans, u_trans);
      const int32_t e_trans_enum = edge_to_int(e_trans, edges);
      edge_transformations[i].push_back(e_trans_enum);
    }
  }
  return edge_transformations;
}

std::vector<edge_inv_t> compute_edge_inversions(
    const std::vector<edge_t>& edges,
    const std::vector<vertex_inv_t>& vertex_inversions, int32_t n) {
  std::vector<edge_inv_t> edge_inversions(vertex_inversions.size());
  for (std::size_t i = 0; i < edge_inversions.size(); ++i) {
    edge_inversions[i].reserve(num_edges(n));
    for (int32_t e = 0; e < num_edges(n); ++e) {
      const vertex_t u_inv = vertex_inversions[i][edges[e].first];
      const vertex_t v_inv = vertex_inversions[i][edges[e].second];
      const edge_t e_inv =
          (u_inv < v_inv) ? edge_t(u_inv, v_inv) : edge_t(v_inv, u_inv);
      const int32_t e_inv_enum = edge_to_int(e_inv, edges);
      edge_inversions[i].push_back(e_inv_enum);
    }
  }
  return edge_inversions;
}

#ifndef N_CUBE_SYMMETRY_H_
#define N_CUBE_SYMMETRY_H_

#include <algorithm>
#include <cstdint>
#include <vector>

#include "common.hpp"
#include "edge.hpp"

/**
 *  Returns all symmetries inherent to the n-cube.
 **/
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

/**
 *  Returns the transformation of a vertex.
 **/
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

/**
 *  Returns the inverse transformation of a vertex.
 **/
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

/**
 *  Returns all inverse vertex transformations.
 *
 *  The element vertex_transformations[i][v] in the returned vector contains the
 *  inverse transformation of vertex v according to the i-th symmetry.
 **/
std::vector<vertex_trans_t> compute_vertex_transformations(int32_t n) {
  const std::vector<symmetry_t> symmetries = compute_symmetries(n);
  std::vector<vertex_trans_t> transformations(symmetries.size());
  for (std::size_t i = 0; i < symmetries.size(); ++i) {
    transformations[i].reserve(num_vertices(n));
    for (vertex_t v = 0; v < num_vertices(n); ++v) {
      transformations[i].push_back(transform_vertex_inv(symmetries[i], v, n));
    }
  }
  return transformations;
}

/**
 *  Returns all inverse edge transformations.
 *
 *  The element edge_transformations[i][e] in the returned vector contains the
 *  inverse transformation of edge e according to the i-th symmetry.
 **/
std::vector<edge_trans_t> compute_edge_transformations(
    const std::vector<edge_t>& edges,
    const std::vector<vertex_trans_t>& vertex_transformations, int32_t n) {
  std::vector<edge_trans_t> edge_transformations(vertex_transformations.size());
  for (std::size_t i = 0; i < edge_transformations.size(); ++i) {
    edge_transformations[i].reserve(num_edges(n));
    for (int32_t e = 0; e < num_edges(n); ++e) {
      const vertex_t u = vertex_transformations[i][edges[e].first];
      const vertex_t v = vertex_transformations[i][edges[e].second];
      const edge_t e_inversion = (u < v) ? edge_t(u, v) : edge_t(v, u);
      edge_transformations[i].push_back(edge_to_int(e_inversion, edges));
    }
  }
  return edge_transformations;
}

#endif  // N_CUBE_SYMMETRY_H_

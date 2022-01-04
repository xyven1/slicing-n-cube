#ifndef N_CUBE_EDGE_H_
#define N_CUBE_EDGE_H_

#include <algorithm>
#include <cstdint>
#include <vector>

#include "common.hpp"

/**
 *  Returns the enumeration of an edge over the lexicographic order of all
 *  edges.
 **/
inline int32_t edge_to_int(const edge_t& e, const std::vector<edge_t>& edges) {
  const auto e_it = std::lower_bound(edges.begin(), edges.end(), e);
  return static_cast<int32_t>(e_it - edges.begin());
}

/**
 *  Returns all edges of the n-cube in lexicographic order.
 **/
std::vector<edge_t> compute_edges(int32_t n) {
  std::vector<edge_t> edges;
  edges.reserve(num_edges(n));
  for (vertex_t v = 0; v < num_vertices(n); ++v) {
    for (int32_t i = 0; i < n; ++i) {
      const vertex_t neighbour = v ^ (1 << i);
      if (v < neighbour) {
        edges.emplace_back(v, neighbour);
      }
    }
  }
  std::sort(edges.begin(), edges.end());
  return edges;
}

#endif  // N_CUBE_EDGE_H_

#ifndef EDGE_H
#define EDGE_H

#include <algorithm>
#include <cstdint>
#include <vector>

#include "common.hpp"

inline int32_t edge_to_int(const edge_t& e, const std::vector<edge_t>& edges) {
  const auto e_it = std::lower_bound(edges.begin(), edges.end(), e);
  return static_cast<int32_t>(e_it - edges.begin());
}

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

#endif

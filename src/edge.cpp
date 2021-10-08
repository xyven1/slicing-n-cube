#include "edge.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

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

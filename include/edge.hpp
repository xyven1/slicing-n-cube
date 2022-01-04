#ifndef N_CUBE_EDGE_H_
#define N_CUBE_EDGE_H_

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

#include "vertex.hpp"

// The two vertices are assumed to be ordered, i.e. e.first < e.second
using edge_t = std::pair<vertex_t, vertex_t>;

// n * 2^n
constexpr int32_t num_edges(int32_t n) { return n << (n - 1); }

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
      const vertex_t u = get_neighbour(v, i);
      if (u < v) {
        edges.emplace_back(u, v);
      }
    }
  }
  std::sort(edges.begin(), edges.end());
  return edges;
}

#endif  // N_CUBE_EDGE_H_

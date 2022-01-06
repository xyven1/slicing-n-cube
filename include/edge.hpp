#ifndef N_CUBE_EDGE_H_
#define N_CUBE_EDGE_H_

#include <algorithm>
#include <array>
#include <cstdint>
#include <utility>
#include <vector>

#include "vertex.hpp"

/* The two vertices are assumed to be ordered, i.e. e.first < e.second. */
using edge_t = std::pair<vertex_t, vertex_t>;

/* The number of edges is n * 2^(n - 1). */
constexpr int32_t num_edges(int32_t n) { return n << (n - 1); }

/**
 *  Returns a symmetric transformation of an edge.
 *
 *  For both vertices every coordinate i is permuted to position permutation[i]
 *  and its sign is flipped if the i-th least significant bit of signs is set.
 **/
template <int32_t N>
edge_t transform_edge(edge_t e, const std::array<int32_t, N>& permutation,
                      int32_t signs) {
  const vertex_t u = transform_vertex<N>(e.first, permutation, signs);
  const vertex_t v = transform_vertex<N>(e.second, permutation, signs);
  return (u < v) ? edge_t(u, v) : edge_t(v, u);
}

/**
 *  Returns the enumeration of an edge over the lexicographic order of all
 *  edges.
 **/
inline int32_t edge_to_int(const edge_t& e, const std::vector<edge_t>& edges) {
  const auto e_it = std::lower_bound(edges.begin(), edges.end(), e);
  return static_cast<int32_t>(e_it - edges.begin());
}

/**
 *  Returns all edges in lexicographic order.
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

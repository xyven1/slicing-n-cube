#ifndef N_CUBE_EDGE_H_
#define N_CUBE_EDGE_H_

#include <algorithm>
#include <array>
#include <cstdint>
#include <utility>

#include "vertex.hpp"

namespace ncube {

/* The two vertices are assumed to be ordered, i.e. e.first < e.second. */
using edge_t = std::pair<vertex_t, vertex_t>;

/* The number of edges is n * 2^(n - 1). */
constexpr int32_t num_edges(int32_t n) { return n << (n - 1); }

/* All edges sorted in lexigraphic order to enable efficient enumeration. */
template <int32_t N>
using edge_lexicon_t = std::array<edge_t, num_edges(N)>;

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
template <int32_t N>
int32_t edge_to_int(const edge_t& e, const edge_lexicon_t<N>& edges) {
  const auto e_it = std::lower_bound(edges.begin(), edges.end(), e);
  return static_cast<int32_t>(e_it - edges.begin());
}

/**
 *  Returns all edges in lexicographic order.
 **/
template <int32_t N>
edge_lexicon_t<N> compute_edges() {
  edge_lexicon_t<N> edges;
  auto it = edges.begin();
  for (vertex_t v = 0; v < num_vertices(N); ++v) {
    for (int32_t i = 0; i < N; ++i) {
      const vertex_t u = get_neighbour(v, i);
      if (u < v) {
        *it++ = edge_t(u, v);
      }
    }
  }
  std::sort(edges.begin(), edges.end());
  return edges;
}

}  // namespace ncube

#endif  // N_CUBE_EDGE_H_

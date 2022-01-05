#ifndef N_CUBE_SYMMETRY_H_
#define N_CUBE_SYMMETRY_H_

#include <array>
#include <cstdint>

#include "edge.hpp"
#include "vertex.hpp"

// n! * 2^n
constexpr int32_t num_symmetries(int32_t n) {
  int32_t factorial = 1;
  for (int32_t i = 2; i <= n; ++i) {
    factorial *= i;
  }
  return (1 << n) * factorial;
}

template <int32_t N>
vertex_t transform_vertex(vertex_t v, const std::array<int32_t, N>& positions,
                          int32_t signs) {
  vertex_t v_trans = 0;
  for (int32_t i = 0; i < N; ++i) {
    const int32_t change_sign = (signs >> i) & 1;
    const int32_t coordinate_bit = (v >> i) & 1;
    const int32_t new_coordinate_bit = coordinate_bit ^ change_sign;
    v_trans |= new_coordinate_bit << positions[i];
  }
  return v_trans;
}

template <int32_t N>
edge_t transform_edge(edge_t e, const std::array<int32_t, N>& positions,
                      int32_t signs) {
  const vertex_t u = transform_vertex<N>(e.first, positions, signs);
  const vertex_t v = transform_vertex<N>(e.second, positions, signs);
  return (u < v) ? edge_t(u, v) : edge_t(v, u);
}

#endif  // N_CUBE_SYMMETRY_H_

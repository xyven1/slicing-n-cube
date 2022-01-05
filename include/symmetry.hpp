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
vertex_t transform_vertex_inv(vertex_t v,
                              const std::array<int32_t, N>& positions,
                              int32_t signs) {
  vertex_t transformation = 0;
  for (int32_t i = 0; i < N; ++i) {
    const int32_t sign_i = (signs >> i) & 1;
    const int32_t val_i = (v >> positions[i]) & 1;
    const int32_t val_i_sign = val_i ^ sign_i;
    const int32_t val_i_sign_and_position = val_i_sign << i;
    transformation |= val_i_sign_and_position;
  }
  return transformation;
}

#endif  // N_CUBE_SYMMETRY_H_

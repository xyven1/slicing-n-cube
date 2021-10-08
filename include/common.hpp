#ifndef COMMON_H
#define COMMON_H

#include <bitset>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

#include "prettyprint.hpp"

template <typename T>
constexpr T pow(T base, T exponent) {
  T power = 1;
  for (T i = 0; i < exponent; ++i) {
    power *= base;
  }
  return power;
}

template <typename T>
constexpr T factorial(T n) {
  T f = 1;
  for (T i = 2; i <= n; ++i) {
    f *= i;
  }
  return f;
}

constexpr int32_t num_vertices(int32_t n) {
    return pow(2, n);
}

constexpr int32_t num_egdes(int32_t n) {
    return n * pow(2, n - 1);
}

constexpr int32_t num_symmetries(int32_t n) {
    return pow(2, n) * factorial(n);
}

constexpr int32_t N = 6;

// The i-th least significant bit stores the i-th coordinate
using vertex_t = int32_t;
using edge_t = std::pair<vertex_t, vertex_t>;
// The MSB stores the sign and the remaining bits represent the new position
using symmetry_t = std::vector<uint32_t>;
// inversion_t[v] stores the inversion of vertex v
using inversion_t = std::vector<vertex_t>;
// transformation_t[v] stores the transformation of vertex v
using transformation_t = std::vector<vertex_t>;
using complex_t = std::bitset<num_vertices(N)>;
using sliceable_set_t = std::bitset<num_egdes(N)>;

#endif

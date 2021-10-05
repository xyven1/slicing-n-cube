#ifndef COMMON_H
#define COMMON_H

#include <bitset>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

#include "prettyprint.hpp"

// The i-th least significant bit stores the i-th coordinate
using vertex_t = uint32_t;
using edge_t = std::pair<vertex_t, vertex_t>;
// The MSB stores the sign and the remaining bits represent the new position
using symmetry_t = std::vector<uint32_t>;
// inversion_t[v] stores the inversion of vertex v
using inversion_t = std::vector<vertex_t>;
// transformation_t[v] stores the transformation of vertex v
using transformation_t = std::vector<vertex_t>;
using complex_t = std::bitset<64>;
using sliceable_set_t = std::bitset<192>;

template <std::size_t N>
bool operator<(const std::bitset<N>& x, const std::bitset<N>& y) {
  for (int i = N - 1; i >= 0; --i) {
    if (x[i] != y[i]) {
      return y[i];
    }
  }
  return false;
}

template <std::size_t N>
bool operator>(const std::bitset<N>& x, const std::bitset<N>& y) {
  for (int i = N - 1; i >= 0; --i) {
    if (x[i] != y[i]) {
      return x[i];
    }
  }
  return false;
}

template <typename T>
T pow(T base, T exponent) {
  T power = 1;
  for (T i = 0; i < exponent; ++i) {
    power *= base;
  }
  return power;
}

template <typename T>
T factorial(T n) {
  T f = 1;
  for (T i = 2; i <= n; ++i) {
    f *= i;
  }
  return f;
}

#endif

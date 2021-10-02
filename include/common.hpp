#ifndef COMMON_H
#define COMMON_H

#include <cstdint>
#include <set>
#include <utility>
#include <vector>

// The i-th least significant bit stores the i-th coordinate
using vertex_t = uint32_t;
using edge_t = std::pair<vertex_t, vertex_t>;
using symmetry_t = std::vector<std::pair<int32_t, int32_t>>;
using complex_t = std::set<vertex_t>;

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

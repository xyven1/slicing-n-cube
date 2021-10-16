#ifndef COMMON_H
#define COMMON_H

#include <bitset>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

#include "prettyprint.hpp"

constexpr int32_t num_vertices(int32_t n) { return 1 << n; }

constexpr int32_t num_edges(int32_t n) { return n << (n - 1); }

constexpr int32_t num_symmetries(int32_t n) {
  int32_t factorial = 1;
  for (int32_t i = 2; i <= n; ++i) {
    factorial *= i;
  }
  return (1 << n) * factorial;
}

template <std::size_t N>
bool operator<(const std::bitset<N>& x, const std::bitset<N>& y) {
  for (std::size_t i = N - 1; i < N; --i) {
    if (x[i] != y[i]) {
      return y[i];
    }
  }
  return false;
}

template <std::size_t N>
bool operator>(const std::bitset<N>& x, const std::bitset<N>& y) {
  for (std::size_t i = N - 1; i < N; --i) {
    if (x[i] != y[i]) {
      return x[i];
    }
  }
  return false;
}

constexpr int32_t N = 5;

// The i-th least significant bit stores the i-th coordinate
using vertex_t = int32_t;

// The two vertices are assumed to be ordered, i.e. edge_t.first < edge_t.second
using edge_t = std::pair<vertex_t, vertex_t>;

// The MSB stores the sign and the remaining bits represent the new position
using symmetry_t = std::vector<uint32_t>;

// vertex_trans_t[v] stores the vertex whose transformation yields v
using vertex_trans_t = std::vector<vertex_t>;

// edge_trans_t[e] stores the edge whose transformation yields e
using edge_trans_t = std::vector<int32_t>;

using complex_t = std::bitset<num_vertices(N)>;

using sliceable_set_t = std::bitset<num_edges(N)>;

#endif

#ifndef EDGE_H
#define EDGE_H

#include <algorithm>
#include <cstdint>
#include <vector>

#include "common.hpp"

inline int32_t edge_to_int(const edge_t& e, const std::vector<edge_t>& edges) {
  const auto e_it = std::lower_bound(edges.begin(), edges.end(), e);
  return static_cast<int32_t>(e_it - edges.begin());
}

std::vector<edge_t> compute_edges(int32_t);

#endif

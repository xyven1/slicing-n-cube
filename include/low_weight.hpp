#ifndef LOW_WEIGHT_H
#define LOW_WEIGHT_H

#include <algorithm>
#include <cstdint>
#include <unordered_set>
#include <vector>

#include "common.hpp"
#include "edge.hpp"
#include "sliceable_set.hpp"

template <int32_t N>
sliceable_set_t<N> low_weight_halfspace_to_sliceable_set(
    const std::vector<int32_t>& normal, double distance,
    const std::vector<edge_t>& edges) {
  sliceable_set_t<N> ss;
  for (const auto& e : edges) {
    int32_t u_scalar = 0, v_scalar = 0;
    for (int32_t i = 0; i < N; ++i) {
      const int32_t u_i = (e.first & (1 << i)) ? 1 : -1;
      const int32_t v_i = (e.second & (1 << i)) ? 1 : -1;
      u_scalar += u_i * normal[i];
      v_scalar += v_i * normal[i];
    }
    if ((u_scalar < distance && v_scalar > distance) ||
        (u_scalar > distance && v_scalar < distance)) {
      ss[edge_to_int(e, edges)] = true;
    }
  }
  return ss;
}

bool next_low_weight_vector(std::vector<int32_t>& halfspace) {
  for (auto it = halfspace.rbegin(); it != halfspace.rend(); ++it) {
    if (*it == -1) {
      *it = 1;
      return true;
    } else {
      *it = -1;
    }
  }
  return false;
}

template <int32_t N>
std::vector<sliceable_set_t<N>> compute_low_weight_sliceable_sets(
    const std::vector<double>& distances, const std::vector<edge_t>& edges) {
  std::unordered_set<sliceable_set_t<N>> sets;
  std::vector<int32_t> normal(N, -1);
  do {
    for (const auto& distance : distances) {
      const sliceable_set_t<N> ss =
          low_weight_halfspace_to_sliceable_set<N>(normal, distance, edges);
      if (ss.any()) {
        sets.insert(ss);
      }
    }
  } while (next_low_weight_vector(normal));
  return std::vector<sliceable_set_t<N>>(sets.begin(), sets.end());
}

template <int32_t N>
int32_t combine_low_weight_sliceable_sets(
    const std::vector<sliceable_set_t<N>>& sets,
    const std::vector<edge_trans_t>& transformations) {
  std::vector<sliceable_set_t<N>> sets_k;
  for (const auto& ss : sets) {
    const auto usr = unique_sliceable_set<N>(ss, transformations);
    if (std::find(sets_k.begin(), sets_k.end(), usr) == sets_k.end()) {
      sets_k.push_back(usr);
    }
  }
  std::cout << "expanded size = " << sets.size() << std::endl;
  std::cout << "k = " << 1 << " size = " << sets_k.size() << std::endl;
  for (int i = 2;; ++i) {
    sets_k = combine_usr_mss<N>(sets_k, sets, transformations);
    std::cout << "k = " << i << " size = " << sets_k.size() << std::endl;
    for (const auto& ss : sets_k) {
      if (ss.all()) {
        return i;
      }
    }
  }
}

#endif

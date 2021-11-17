#include "low_weight.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

#include "common.hpp"
#include "edge.hpp"
#include "sliceable_set.hpp"

template <int32_t N>
sliceable_set_t<N> low_weight_halfspace_to_sliceable_set(
    const std::vector<int32_t>& normal, int32_t distance,
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

bool next_low_weight_vector(std::vector<int32_t>& halfspace, int32_t min_weight,
                            int32_t max_weight) {
  for (auto it = halfspace.rbegin(); it != halfspace.rend(); ++it) {
    if (*it == max_weight) {
      *it = min_weight;
    } else if (*it == -1) {
      *it = 1;
      return true;
    } else {
      *it += 1;
      return true;
    }
  }
  return false;
}

template <int32_t N>
std::vector<sliceable_set_t<N>> compute_low_weight_sliceable_sets(
    int32_t min_weight, int32_t max_weight, const std::vector<edge_t>& edges) {
  std::vector<sliceable_set_t<N>> sets;
  std::vector<int32_t> normal(N, min_weight);
  do {
    for (int i = min_weight; i <= max_weight; ++i) {
      const sliceable_set_t<N> ss =
          low_weight_halfspace_to_sliceable_set<N>(normal, i, edges);
      if (ss.any()) {
        // std::cout << normal << " " << i << " " << ss << std::endl;
        sets.push_back(ss);
      }
    }
  } while (next_low_weight_vector(normal, min_weight, max_weight));
  std::sort(sets.begin(), sets.end());
  sets.erase(std::unique(sets.begin(), sets.end()), sets.end());
  sets.shrink_to_fit();
  return sets;
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

template sliceable_set_t<2> low_weight_halfspace_to_sliceable_set<2>(
    const std::vector<int32_t>& normal, int32_t distance,
    const std::vector<edge_t>& edges);
template sliceable_set_t<3> low_weight_halfspace_to_sliceable_set<3>(
    const std::vector<int32_t>& normal, int32_t distance,
    const std::vector<edge_t>& edges);
template sliceable_set_t<4> low_weight_halfspace_to_sliceable_set<4>(
    const std::vector<int32_t>& normal, int32_t distance,
    const std::vector<edge_t>& edges);
template sliceable_set_t<5> low_weight_halfspace_to_sliceable_set<5>(
    const std::vector<int32_t>& normal, int32_t distance,
    const std::vector<edge_t>& edges);
template sliceable_set_t<6> low_weight_halfspace_to_sliceable_set<6>(
    const std::vector<int32_t>& normal, int32_t distance,
    const std::vector<edge_t>& edges);
template sliceable_set_t<7> low_weight_halfspace_to_sliceable_set<7>(
    const std::vector<int32_t>& normal, int32_t distance,
    const std::vector<edge_t>& edges);

template std::vector<sliceable_set_t<2>> compute_low_weight_sliceable_sets<2>(
    int32_t min_weight, int32_t max_weight, const std::vector<edge_t>& edges);
template std::vector<sliceable_set_t<3>> compute_low_weight_sliceable_sets<3>(
    int32_t min_weight, int32_t max_weight, const std::vector<edge_t>& edges);
template std::vector<sliceable_set_t<4>> compute_low_weight_sliceable_sets<4>(
    int32_t min_weight, int32_t max_weight, const std::vector<edge_t>& edges);
template std::vector<sliceable_set_t<5>> compute_low_weight_sliceable_sets<5>(
    int32_t min_weight, int32_t max_weight, const std::vector<edge_t>& edges);
template std::vector<sliceable_set_t<6>> compute_low_weight_sliceable_sets<6>(
    int32_t min_weight, int32_t max_weight, const std::vector<edge_t>& edges);
template std::vector<sliceable_set_t<7>> compute_low_weight_sliceable_sets<7>(
    int32_t min_weight, int32_t max_weight, const std::vector<edge_t>& edges);

template int32_t combine_low_weight_sliceable_sets<2>(
    const std::vector<sliceable_set_t<2>>& sets,
    const std::vector<edge_trans_t>& transformations);
template int32_t combine_low_weight_sliceable_sets<3>(
    const std::vector<sliceable_set_t<3>>& sets,
    const std::vector<edge_trans_t>& transformations);
template int32_t combine_low_weight_sliceable_sets<4>(
    const std::vector<sliceable_set_t<4>>& sets,
    const std::vector<edge_trans_t>& transformations);
template int32_t combine_low_weight_sliceable_sets<5>(
    const std::vector<sliceable_set_t<5>>& sets,
    const std::vector<edge_trans_t>& transformations);
template int32_t combine_low_weight_sliceable_sets<6>(
    const std::vector<sliceable_set_t<6>>& sets,
    const std::vector<edge_trans_t>& transformations);
template int32_t combine_low_weight_sliceable_sets<7>(
    const std::vector<sliceable_set_t<7>>& sets,
    const std::vector<edge_trans_t>& transformations);

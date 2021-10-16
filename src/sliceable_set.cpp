#include "sliceable_set.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

#include "common.hpp"

bool is_subset(const sliceable_set_t& subset,
               const std::vector<sliceable_set_t>& combos) {
  for (const sliceable_set_t& combo : combos) {
    if ((subset | combo) == combo) {
      return true;
    }
  }
  return false;
}

sliceable_set_t unique_sliceable_set(
    const sliceable_set_t& ss, const std::vector<edge_trans_t>& transformations,
    int32_t n) {
  sliceable_set_t min_ss(ss);
  for (const edge_trans_t& transformation : transformations) {
    sliceable_set_t ss_trans;
    bool is_new_min = false;
    for (int32_t e = num_edges(n) - 1; e >= 0; --e) {
      const int32_t e_inversion = transformation[e];
      ss_trans[e] = ss[e_inversion];
      if (!is_new_min && ss_trans[e] && !min_ss[e]) {
        // min_ss is still smaller
        break;
      }
      is_new_min |= !ss_trans[e] && min_ss[e];
    }
    if (is_new_min) {
      min_ss = ss_trans;
    }
  }
  return min_ss;
}

std::vector<sliceable_set_t> combine_usr_mss(
    const std::vector<sliceable_set_t>& usr,
    const std::vector<sliceable_set_t>& mss,
    const std::vector<edge_trans_t>& transformations, int32_t n) {
  std::vector<sliceable_set_t> combos;
  for (const sliceable_set_t& set_1 : usr) {
    for (const sliceable_set_t& set_2 : mss) {
      sliceable_set_t combo = set_1 | set_2;
      combo = unique_sliceable_set(combo, transformations, n);
      if (!is_subset(combo, combos)) {
        const auto p = [combo](const sliceable_set_t& ss) {
          return (ss | combo) == combo;
        };
        const auto it = remove_if(combos.begin(), combos.end(), p);
        combos.erase(it, combos.end());
        combos.push_back(combo);
      }
    }
  }
  return combos;
}

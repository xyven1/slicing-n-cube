#ifndef N_CUBE_SLICE_CUBE_H_
#define N_CUBE_SLICE_CUBE_H_

#include <algorithm>
#include <cstdint>
#include <vector>

#include "edge.hpp"
#include "sliceable_set.hpp"

namespace ncube {

/**
 *  Returns the smallest k so that a combination of k sliceable sets given by
 *  their unique symmetric representations slices the n-cube.
 *
 *  If k is larger than 7, -1 is returned instead.
 */
template <int32_t N>
int32_t slice_cube_min(const std::vector<sliceable_set_t<N>>& usr_1,
                       const edge_lexicon_t<N>& edges) {
  const auto all = [](const sliceable_set_t<N>& ss) { return ss.all(); };
  if (std::any_of(usr_1.begin(), usr_1.end(), all)) {
    return 1;
  }
  const auto mss_1 = expand_usr<N>(usr_1, edges);
  if (pairwise_unions_slice_cube<N>(usr_1, mss_1)) {
    return 2;
  }
  const auto usr_2 = pairwise_unions<N>(usr_1, mss_1, edges);
  if (pairwise_unions_slice_cube<N>(usr_2, mss_1)) {
    return 3;
  }
  const auto mss_2 = expand_usr<N>(usr_2, edges);
  if (pairwise_unions_slice_cube<N>(usr_2, mss_2)) {
    return 4;
  }
  const auto usr_3 = pairwise_unions<N>(usr_2, mss_1, edges);
  if (pairwise_unions_slice_cube<N>(usr_3, mss_2)) {
    return 5;
  }
  const auto mss_3 = expand_usr<N>(usr_3, edges);
  if (pairwise_unions_slice_cube<N>(usr_3, mss_3)) {
    return 6;
  }
  const auto usr_4 = pairwise_unions<N>(usr_2, mss_2, edges);
  if (pairwise_unions_slice_cube<N>(usr_4, mss_3)) {
    return 7;
  }
  return -1;
}

}  // namespace ncube

#endif  // N_CUBE_SLICE_CUBE_H_

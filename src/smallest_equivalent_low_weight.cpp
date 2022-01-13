#include <cstdint>
#include <iostream>

#include "complex.hpp"
#include "edge.hpp"
#include "low_weight.hpp"
#include "slice_cube.hpp"
#include "sliceable_set.hpp"
#include "vertex.hpp"

using namespace ncube;

/**
 *  Returns the smallest i so that low weight halfspaces whose normal vector
 *  contains only values in {-i, ..., i} have the same maximal sliceable sets as
 *  in the general setting.
 **/
template <int32_t N>
int32_t equivalent_low_weight_mss() {
  const auto complexes = compute_complexes<N>(is_complex_degree_one<N>);
  const auto edges = compute_edges<N>();
  const auto usr = complexes_to_usr<N>(complexes, edges);
  const auto mss = expand_usr<N>(usr, edges);
  for (int i = 0;; ++i) {
    auto mss_low_weight = compute_low_weight_mss<N>(i, edges);
    std::sort(mss_low_weight.begin(), mss_low_weight.end());
    if (mss_low_weight == mss) {
      return i;
    }
  }
}

/**
 *  Returns the smallest i so that low weight halfspaces whose normal vector
 *  contains only values in {-i, ..., i} slice the n-cube with the same number
 *  of sliceable sets as in the general setting.
 **/
template <int32_t N>
int32_t equivalent_low_weight_slice_cube_min() {
  const auto complexes = compute_complexes<N>(is_complex_degree_one<N>);
  const auto edges = compute_edges<N>();
  const auto usr = complexes_to_usr<N>(complexes, edges);
  const auto k = slice_cube_min<N>(usr, edges);
  for (int i = 0;; ++i) {
    const auto mss_low_weight = compute_low_weight_mss<N>(i, edges);
    const auto usr_low_weight = reduce_to_usr<N>(mss_low_weight, edges);
    const auto k_low_weight = slice_cube_min<N>(usr_low_weight, edges);
    if (k == k_low_weight) {
      return i;
    }
  }
}

int main() {
  std::cout << "n=3 " << equivalent_low_weight_mss<3>() << std::endl;
  std::cout << "n=4 " << equivalent_low_weight_mss<4>() << std::endl;
  std::cout << "n=5 " << equivalent_low_weight_mss<5>() << std::endl;
  std::cout << "n=3 " << equivalent_low_weight_slice_cube_min<3>() << std::endl;
  std::cout << "n=4 " << equivalent_low_weight_slice_cube_min<4>() << std::endl;
  std::cout << "n=5 " << equivalent_low_weight_slice_cube_min<5>() << std::endl;
}

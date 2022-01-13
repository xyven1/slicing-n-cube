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
 *  Outputs the smallest i so that low weight halfspaces whose normal vector
 *  contains only values in {-i, ..., i} have the same maximal sliceable sets as
 *  in the general setting.
 **/
template <int32_t N>
void equivalent_low_weight_mss() {
  std::cout << "n = " << N << std::endl;
  const auto complexes = compute_complexes<N>(is_complex_degree_one<N>);
  const auto edges = compute_edges<N>();
  const auto usr = complexes_to_usr<N>(complexes, edges);
  const auto mss = expand_usr<N>(usr, edges);
  std::cout << "  |mss| = " << mss.size() << std::endl;
  for (int i = 1;; ++i) {
    const auto mss_low_weight = compute_low_weight_mss<N>(i, edges);
    std::cout << "  |mss_" << i << "| = " << mss_low_weight.size() << std::endl;
    if (mss_low_weight == mss) {
      std::cout << "  smallest i to have equivalent mss is " << i << std::endl;
      return;
    }
  }
}

/**
 *  Outputs the smallest i so that low weight halfspaces whose normal vector
 *  contains only values in {-i, ..., i} slice the n-cube with the same number
 *  of sliceable sets as in the general setting.
 **/
template <int32_t N>
void equivalent_low_weight_slice_cube_min() {
  std::cout << "n = " << N << std::endl;
  const auto complexes = compute_complexes<N>(is_complex_degree_one<N>);
  const auto edges = compute_edges<N>();
  const auto usr = complexes_to_usr<N>(complexes, edges);
  const auto k = slice_cube_min<N>(usr, edges);
  std::cout << "  k = " << k << std::endl;
  for (int i = 1;; ++i) {
    const auto mss_low_weight = compute_low_weight_mss<N>(i, edges);
    const auto usr_low_weight = reduce_to_usr<N>(mss_low_weight, edges);
    const auto k_low_weight = slice_cube_min<N>(usr_low_weight, k, edges);
    std::cout << "  k_" << i << " = " << k << std::endl;
    if (k == k_low_weight) {
      std::cout << "  smallest i to slice the n-cube with the same number of "
                   "hyperplanes is "
                << i << std::endl;
      return;
    }
  }
}

int main() {
  equivalent_low_weight_mss<2>();
  equivalent_low_weight_mss<3>();
  equivalent_low_weight_mss<4>();
  equivalent_low_weight_mss<5>();
  equivalent_low_weight_slice_cube_min<2>();
  equivalent_low_weight_slice_cube_min<3>();
  equivalent_low_weight_slice_cube_min<4>();
  equivalent_low_weight_slice_cube_min<5>();
}

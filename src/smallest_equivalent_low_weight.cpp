#include <cstdint>
#include <iostream>

#include "complex.hpp"
#include "edge.hpp"
#include "low_weight.hpp"
#include "sliceable_set.hpp"
#include "vertex.hpp"

using namespace ncube;

/**
 *  Returns the smallest k so that low weight halfspaces whose normal vector
 *  contains only values in {-k, ..., k} cover all maximal sliceable sets.
 **/
template <int32_t N>
int32_t smallest_equivalent_low_weight() {
  const auto complexes = compute_complexes<N>(is_complex_degree_one<N>);
  const auto edges = compute_edges<N>();
  const auto usr = complexes_to_usr<N>(complexes, edges);
  const auto mss = expand_usr<N>(usr, edges);
  for (int k = 0;; ++k) {
    const auto usr_low_weight = compute_low_weight_usr_mss<N>(k, edges);
    const auto mss_low_weight = expand_usr<N>(usr_low_weight, edges);
    if (mss_low_weight == mss) {
      return k;
    }
    std::cout << "k = " << k << ", |mss| = " << mss.size() << "|mss_low_weight|"
              << mss_low_weight.size() << std::endl;
  }
}

int main() {
  std::cout << "n = 3\n" << smallest_equivalent_low_weight<3>() << std::endl;
  std::cout << "n = 4\n" << smallest_equivalent_low_weight<4>() << std::endl;
  std::cout << "n = 5\n" << smallest_equivalent_low_weight<5>() << std::endl;
  std::cout << "n = 6\n" << smallest_equivalent_low_weight<6>() << std::endl;
}

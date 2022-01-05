#include <cstdint>
#include <iostream>

#include "complex.hpp"
#include "edge.hpp"
#include "low_weight.hpp"
#include "sliceable_set.hpp"
#include "symmetry.hpp"
#include "vertex.hpp"

template <int32_t N>
int32_t smallest_low_weight() {
  const auto complexes = compute_cut_complexes_degree_one<N>();
  const auto edges = compute_edges(N);
  const auto usr = complexes_to_usr<N>(complexes, edges);
  const auto mss = usr_to_mss<N>(usr, edges);
  for (int k = 0;; ++k) {
    auto mss_low_weight = compute_low_weight_mss<N>(edges, k);
    std::sort(mss_low_weight.begin(), mss_low_weight.end());
    if (mss_low_weight == mss) {
      return k;
    }
    std::cout << "k = " << k << ", |mss| = " << mss.size() << "|mss_low_weight|"
              << mss_low_weight.size() << std::endl;
  }
}

int main() {
  std::cout << "n = 3\n" << smallest_low_weight<3>() << std::endl;
  std::cout << "n = 4\n" << smallest_low_weight<4>() << std::endl;
  std::cout << "n = 5\n" << smallest_low_weight<5>() << std::endl;
  std::cout << "n = 6\n" << smallest_low_weight<6>() << std::endl;
}

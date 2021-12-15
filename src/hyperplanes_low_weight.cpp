#include <cstdint>
#include <iostream>

#include "common.hpp"
#include "complex.hpp"
#include "edge.hpp"
#include "low_weight.hpp"
#include "symmetry.hpp"

template <int32_t N>
auto compute_sorted_mss() {
  const auto vertex_transformations = compute_vertex_transformations(N);
  const auto complexes = compute_cut_complexes<N>(vertex_transformations);
  const auto edges = compute_edges(N);
  auto mss = complexes_to_mss<N>(complexes, vertex_transformations, edges);
  std::sort(mss.begin(), mss.end());
  return mss;
}

template <int32_t N>
int32_t smallest_low_weight() {
  const auto mss = compute_sorted_mss<N>();
  const auto edges = compute_edges(N);
  for (int k = 0;; ++k) {
    auto mss_low_weight = compute_low_weight_mss<N>(edges, k);
    std::sort(mss_low_weight.begin(), mss_low_weight.end());
    if (mss_low_weight == mss) {
      return k;
    }
    std::cout << "k = " << k << ", |mss| = " << mss.size() << "|mss_low_weight|" << mss_low_weight.size() << std::endl;
  }
}

int main() {
  std::cout << smallest_low_weight<3>() << std::endl;
  std::cout << smallest_low_weight<4>() << std::endl;
  std::cout << smallest_low_weight<5>() << std::endl;
  std::cout << smallest_low_weight<6>() << std::endl;
}

#include <cstdint>
#include <iostream>

#include "complex.hpp"
#include "edge.hpp"
#include "sliceable_set.hpp"
#include "vertex.hpp"

using namespace ncube;

void slice_cube_5() {
  constexpr int32_t N = 5;
  const auto edges = compute_edges<N>();
  const auto complexes = compute_complexes<N>(is_complex_degree_one<N>);
  const auto usr_1 = complexes_to_usr<N>(complexes, edges);
  std::cout << usr_1.size() << std::endl;
  const auto mss_1 = expand_usr<N>(usr_1, edges);
  std::cout << mss_1.size() << std::endl;
  const auto usr_2 = pairwise_unions<N>(usr_1, mss_1, edges);
  std::cout << usr_2.size() << std::endl;
  const auto mss_2 = expand_usr<N>(usr_2, edges);
  std::cout << mss_2.size() << std::endl;
  const auto slices_cube = pairwise_unions_slice_cube<N>(usr_2, mss_2);
  std::cout << slices_cube << std::endl;
}

int main() { slice_cube_5(); }

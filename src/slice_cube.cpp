#include <cstdint>
#include <iostream>

#include "complex.hpp"
#include "edge.hpp"
#include "sliceable_set.hpp"
#include "symmetry.hpp"
#include "vertex.hpp"

bool slice_cube_5() {
  constexpr int32_t N = 5;
  const auto edges = compute_edges(N);
  const auto complexes = compute_cut_complexes_degree_one<N>();
  const auto usr_1 = complexes_to_usr<N>(complexes, edges);
  const auto mss_1 = usr_to_mss<N>(usr_1, edges);
  const auto usr_2 = pairwise_unions<N>(usr_1, mss_1, edges);
  const auto mss_2 = usr_to_mss<N>(usr_2, edges);
  return pairwise_unions_slice_cube<N>(usr_2, mss_2);
}

int main() { std::cout << slice_cube_5() << std::endl; }

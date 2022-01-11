#include <cstdint>
#include <iostream>

#include "complex.hpp"
#include "edge.hpp"
#include "slice_cube.hpp"
#include "sliceable_set.hpp"
#include "vertex.hpp"

using namespace ncube;

template <int32_t N>
int32_t slice_cube_min_degree_two() {
  const auto edges = compute_edges<N>();
  const auto complexes = compute_complexes<N>(is_complex_degree_two<N>);
  const auto usr = complexes_to_usr<N>(complexes, edges);
  const auto k = slice_cube_min<N>(usr, edges);
  return k;
}

bool slice_5_cube_with_2_hyperplanes() {
  constexpr int32_t N = 5;
  const auto edges = compute_edges<N>();
  const auto complexes = compute_complexes<N>(is_complex_degree_two<N>);
  const auto usr = complexes_to_usr<N>(complexes, edges);
  const auto mss = expand_usr<N>(usr, edges);
  const auto slices_cube = pairwise_unions_slice_cube<N>(usr, mss);
  return slices_cube;
}

bool slice_5_cube_with_3_hyperplanes() {
  constexpr int32_t N = 5;
  const auto edges = compute_edges<N>();
  const auto complexes = compute_complexes<N>(is_complex_degree_two<N>);
  const auto usr = complexes_to_usr<N>(complexes, edges);
  const auto mss = expand_usr<N>(usr, edges);
  for (auto it_1 = mss.rbegin(); it_1 != mss.rend(); ++it_1) {
    for (auto it_2 = mss.rbegin(); it_2 != mss.rend(); ++it_2) {
      for (auto it_3 = usr.rbegin(); it_3 != usr.rend(); ++it_3) {
        if ((*it_1 | *it_2 | *it_3).all()) {
          return true;
        }
      }
    }
  }
  return false;
}

int main() {
  std::cout << "Minimum number of degree two polynomials to slice the 2-cube: "
            << slice_cube_min_degree_two<2>() << std::endl;
  std::cout << "Minimum number of degree two polynomials to slice the 3-cube: "
            << slice_cube_min_degree_two<3>() << std::endl;
  std::cout << "Minimum number of degree two polynomials to slice the 4-cube: "
            << slice_cube_min_degree_two<4>() << std::endl;
  std::cout << "Can two degree two polynomials slice the 5-cube: "
            << slice_5_cube_with_2_hyperplanes() << std::endl;
  std::cout << "Can three degree two polynomials slice the 5-cube: "
            << slice_5_cube_with_3_hyperplanes() << std::endl;
}

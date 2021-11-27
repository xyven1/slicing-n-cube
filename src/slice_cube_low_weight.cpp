#include <cstdint>
#include <iostream>

#include "common.hpp"
#include "edge.hpp"
#include "low_weight.hpp"
#include "symmetry.hpp"

template <int32_t N>
int32_t slice_cube_low_weight(const std::vector<double>& distances) {
  const auto edges = compute_edges(N);
  const auto vertex_transformations = compute_vertex_transformations(N);
  const auto edge_transformations =
      compute_edge_transformations(edges, vertex_transformations, N);
  const auto low_weight_sets =
      compute_low_weight_sliceable_sets<N>(distances, edges);
  const auto k = combine_low_weight_sliceable_sets<N>(low_weight_sets,
                                                      edge_transformations);
  return k;
}

int main() {
  constexpr int32_t N = 5;
  std::vector<double> distances = {-1, 0, 1};
  // for (double i = -N; i <= N; i += 0.5) {
  //   distances.push_back(i);
  // }
  std::cout << slice_cube_low_weight<N>(distances) << std::endl;
}

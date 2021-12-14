#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_set>

#include "common.hpp"
#include "complex.hpp"
#include "edge.hpp"
#include "low_weight.hpp"
#include "sliceable_set.hpp"
#include "symmetry.hpp"

int main() {
  constexpr int32_t N = 5;
  const auto edges = compute_edges(N);
  const auto vertex_transformations = compute_vertex_transformations(N);
  const auto edge_transformations =
      compute_edge_transformations(edges, vertex_transformations, N);
  std::vector<double> distances;
  for (double i = -N; i <= N; i += 0.5) {
    distances.push_back(i);
  }
  const auto low_weight_sets =
      compute_one_weight_sliceable_sets<N>(distances, edges);
  const auto k = combine_low_weight_sliceable_sets<N>(low_weight_sets,
                                                      edge_transformations);
  std::cout << k << std::endl;
}

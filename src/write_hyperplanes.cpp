#include <cstdint>
#include <iostream>
#include <vector>

#include "common.hpp"
#include "edge.hpp"
#include "low_weight.hpp"

template <int32_t N>
void write_one_weight_halfspaces() {
  const auto edges = compute_edges(N);
  std::vector<double> distances;
  for (double d = 0.0; d < N; d += 0.5) {
    distances.push_back(d);
  }
  write_one_weight_halfspaces_to_file<N>(
      distances, edges,
      NCUBE_DIR "one_weight_halfspaces_all_" + std::to_string(N) + ".txt");
  write_one_weight_halfspaces_to_file<N>(
      {0, 1}, edges,
      NCUBE_DIR "one_weight_halfspaces_one_" + std::to_string(N) + ".txt");
}

int main() {
  write_one_weight_halfspaces<3>();
  write_one_weight_halfspaces<4>();
  write_one_weight_halfspaces<5>();
  write_one_weight_halfspaces<6>();
  write_one_weight_halfspaces<7>();
}

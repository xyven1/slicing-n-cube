#include <cstdint>
#include <iostream>
#include <vector>

#include "edge.hpp"
#include "low_weight.hpp"
#include "vertex.hpp"

template <int32_t N>
void write_one_weight_halfspaces() {
  const auto edges = compute_edges(N);
  std::vector<int32_t> distances;
  for (int32_t i = 0; i < N; ++i) {
    distances.push_back(i);
  }
  const auto file_name_any = "any_threshold_" + std::to_string(N) + ".txt";
  const auto file_name_one = "one_threshold_" + std::to_string(N) + ".txt";
  const auto file_path_any = NCUBE_DIR "one_weight/" + file_name_any;
  const auto file_path_one = NCUBE_DIR "one_weight/" + file_name_one;
  write_one_weight_halfspaces_to_file<N>(distances, edges, file_path_any);
  write_one_weight_halfspaces_to_file<N>({0, 1}, edges, file_path_one);
}

template <int32_t N>
void write_low_weight_halfspaces(int32_t max) {
  const auto edges = compute_edges(N);
  const auto file_name =
      "max_" + std::to_string(max) + "_" + std::to_string(N) + ".txt";
  const auto file_path = NCUBE_DIR "low_weight/" + file_name;
  write_low_weight_halfspaces_to_file<N>(max, edges, file_path);
}

int main() {
  write_one_weight_halfspaces<3>();
  write_one_weight_halfspaces<4>();
  write_one_weight_halfspaces<5>();
  write_one_weight_halfspaces<6>();
  write_one_weight_halfspaces<7>();
  for (int32_t i = 1; i <= 5; ++i) {
    write_low_weight_halfspaces<3>(i);
    write_low_weight_halfspaces<4>(i);
    write_low_weight_halfspaces<5>(i);
    write_low_weight_halfspaces<6>(i);
    write_low_weight_halfspaces<7>(i);
  }
}

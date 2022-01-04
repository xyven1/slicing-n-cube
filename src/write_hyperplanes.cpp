#include <cstdint>
#include <filesystem>
#include <iostream>
#include <string>
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
  constexpr auto dir = N_CUBE_OUT_DIR "/one_weight";
  std::filesystem::create_directories(dir);
  const auto path_any = dir + ("/any_threshold_" + std::to_string(N) + ".txt");
  const auto path_one = dir + ("/one_threshold_" + std::to_string(N) + ".txt");
  write_one_weight_halfspaces_to_file<N>(distances, edges, path_any);
  write_one_weight_halfspaces_to_file<N>({0, 1}, edges, path_one);
}

template <int32_t N>
void write_low_weight_halfspaces(int32_t max) {
  const auto edges = compute_edges(N);
  constexpr auto dir = N_CUBE_OUT_DIR "/one_weight";
  std::filesystem::create_directories(dir);
  const auto path =
      dir + ("/max_" + std::to_string(max) + "_" + std::to_string(N) + ".txt");
  write_low_weight_halfspaces_to_file<N>(max, edges, path);
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

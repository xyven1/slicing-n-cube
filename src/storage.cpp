#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <string>

#include "complex.hpp"
#include "edge.hpp"
#include "sliceable_set.hpp"
#include "symmetry.hpp"
#include "vertex.hpp"

template <int32_t N>
void write_degree_two_sliceable_sets() {
  const auto vertex_transformations = compute_vertex_transformations(N);
  const auto complexes = compute_cut_complexes_degree_two<N>();
  const auto edges = compute_edges(N);
  const auto edge_transformations =
      compute_edge_transformations(edges, vertex_transformations, N);
  auto usr = complexes_to_usr<N>(complexes, edges);
  auto mss = usr_to_mss<N>(usr, edges);
  constexpr auto dir = N_CUBE_OUT_DIR "/degree_two";
  std::filesystem::create_directories(dir);
  const auto path_usr_1 = dir + ("/" + std::to_string(N) + "_usr_1.bin");
  const auto path_mss_1 = dir + ("/" + std::to_string(N) + "_mss_1.bin");
  write_to_file<N>(usr, path_usr_1);
  write_to_file<N>(mss, path_mss_1);
}

template <int32_t N>
void write_degree_one_sliceable_sets() {
  const auto vertex_transformations = compute_vertex_transformations(N);
  const auto complexes = compute_cut_complexes_degree_one<N>();
  const auto edges = compute_edges(N);
  const auto edge_transformations =
      compute_edge_transformations(edges, vertex_transformations, N);
  auto usr = complexes_to_usr<N>(complexes, edges);
  auto mss = usr_to_mss<N>(usr, edges);
  auto usr_2 = combine_usr_mss<N>(usr, mss, edges);
  auto mss_2 = usr_to_mss<N>(usr_2, edges);
  constexpr auto dir = N_CUBE_OUT_DIR "/degree_one";
  std::filesystem::create_directories(dir);
  const auto path_usr_1 = dir + ("/" + std::to_string(N) + "_usr_1.bin");
  const auto path_mss_1 = dir + ("/" + std::to_string(N) + "_mss_1.bin");
  const auto path_usr_2 = dir + ("/" + std::to_string(N) + "_usr_2.bin");
  const auto path_mss_2 = dir + ("/" + std::to_string(N) + "_mss_2.bin");
  write_to_file<N>(usr, path_usr_1);
  write_to_file<N>(mss, path_mss_1);
  write_to_file<N>(usr_2, path_usr_2);
  write_to_file<N>(mss_2, path_mss_2);
}

int main() {
  write_degree_one_sliceable_sets<2>();
  write_degree_one_sliceable_sets<3>();
  write_degree_one_sliceable_sets<4>();
  write_degree_one_sliceable_sets<5>();
  write_degree_two_sliceable_sets<2>();
  write_degree_two_sliceable_sets<3>();
  write_degree_two_sliceable_sets<4>();
  write_degree_two_sliceable_sets<5>();
}

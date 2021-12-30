#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>

#include "common.hpp"
#include "complex.hpp"
#include "edge.hpp"
#include "sliceable_set.hpp"
#include "symmetry.hpp"

template <int32_t N>
void write_degree_two_sliceable_sets() {
  const auto vertex_transformations = compute_vertex_transformations(N);
  const auto complexes =
      compute_cut_complexes_degree_two<N>(vertex_transformations);
  const auto edges = compute_edges(N);
  const auto edge_transformations =
      compute_edge_transformations(edges, vertex_transformations, N);
  auto usr = complexes_to_usr<N>(complexes, edge_transformations, edges);
  auto mss = complexes_to_mss<N>(complexes, vertex_transformations, edges);
  std::cout << usr.size() << std::endl;
  std::cout << mss.size() << std::endl;
  std::sort(usr.begin(), usr.end());
  std::sort(mss.begin(), mss.end());
  write_to_file<N>(usr,
                   NCUBE_DIR "degree_two/" + std::to_string(N) + "_usr_1.bin");
  write_to_file<N>(mss,
                   NCUBE_DIR "degree_two/" + std::to_string(N) + "_mss_1.bin");
}

template <int32_t N>
void write_two_sliceable_sets() {
  const auto vertex_transformations = compute_vertex_transformations(N);
  const auto complexes =
      compute_cut_complexes_degree_one<N>(vertex_transformations);
  const auto edges = compute_edges(N);
  const auto edge_transformations =
      compute_edge_transformations(edges, vertex_transformations, N);
  auto usr = complexes_to_usr<N>(complexes, edge_transformations, edges);
  auto mss = complexes_to_mss<N>(complexes, vertex_transformations, edges);
  std::cout << usr.size() << std::endl;
  std::cout << mss.size() << std::endl;
  auto usr_2 = combine_usr_mss<N>(usr, mss, edge_transformations);
  std::cout << usr_2.size() << std::endl;
  auto mss_2 = usr_to_mss<N>(usr_2, edge_transformations);
  std::cout << mss_2.size() << std::endl;
  std::sort(usr.begin(), usr.end());
  std::sort(mss.begin(), mss.end());
  std::sort(usr_2.begin(), usr_2.end());
  std::sort(mss_2.begin(), mss_2.end());
  write_to_file<N>(usr, NCUBE_DIR + std::to_string(N) + "_usr_1.bin");
  write_to_file<N>(mss, NCUBE_DIR + std::to_string(N) + "_mss_1.bin");
  write_to_file<N>(usr_2, NCUBE_DIR + std::to_string(N) + "_usr_2.bin");
  write_to_file<N>(mss_2, NCUBE_DIR + std::to_string(N) + "_mss_2.bin");
}

int main() {
  write_two_sliceable_sets<2>();
  write_two_sliceable_sets<3>();
  write_two_sliceable_sets<4>();
  write_two_sliceable_sets<5>();
  write_degree_two_sliceable_sets<2>();
  write_degree_two_sliceable_sets<3>();
  write_degree_two_sliceable_sets<4>();
  write_degree_two_sliceable_sets<5>();
}

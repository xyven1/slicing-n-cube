#include <cstdint>
#include <iostream>
#include <string>

#include "common.hpp"
#include "complex.hpp"
#include "edge.hpp"
#include "sliceable_set.hpp"
#include "symmetry.hpp"

template <int32_t N>
void write_two_sliceable_sets() {
  const std::vector<vertex_trans_t> vertex_transformations =
      compute_vertex_transformations(N);
  const std::vector<complex_t<N>> complexes =
      compute_cut_complexes<N>(vertex_transformations);
  const std::vector<edge_t> edges = compute_edges(N);
  const std::vector<edge_trans_t> edge_transformations =
      compute_edge_transformations(edges, vertex_transformations, N);
  const std::vector<sliceable_set_t<N>> usr =
      complexes_to_usr<N>(complexes, edge_transformations, edges);
  const std::vector<sliceable_set_t<N>> mss =
      complexes_to_mss<N>(complexes, vertex_transformations, edges);
  std::cout << usr.size() << std::endl;
  std::cout << mss.size() << std::endl;
  const std::vector<sliceable_set_t<N>> usr_2 =
      combine_usr_mss<N>(usr, mss, edge_transformations);
  std::cout << usr_2.size() << std::endl;
  const std::vector<sliceable_set_t<N>> mss_2 =
      usr_to_mss<N>(usr_2, edge_transformations);
  std::cout << mss_2.size() << std::endl;
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
}

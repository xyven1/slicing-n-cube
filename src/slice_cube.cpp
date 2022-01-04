#include <cstdint>
#include <iostream>

#include "complex.hpp"
#include "edge.hpp"
#include "sliceable_set.hpp"
#include "symmetry.hpp"
#include "vertex.hpp"

bool slice_cube_5() {
  constexpr int32_t N = 5;
  const auto vertex_transformations = compute_vertex_transformations(N);
  const auto complexes = compute_cut_complexes_degree_one<N>();
  const auto edges = compute_edges(N);
  const auto edge_transformations =
      compute_edge_transformations(edges, vertex_transformations, N);
  const auto usr = complexes_to_usr<N>(complexes, edge_transformations, edges);
  const auto mss = usr_to_mss<N>(usr, edge_transformations);
  const auto usr_2 = combine_usr_mss<N>(usr, mss, edge_transformations);
  const auto mss_2 = usr_to_mss<N>(usr_2, edge_transformations);
  return combine_usr_mss_final<N>(usr, mss);
}

int main() { std::cout << slice_cube_5() << std::endl; }

#include <cstdint>
#include <iostream>

#include "complex.hpp"
#include "edge.hpp"
#include "sliceable_set.hpp"
#include "symmetry.hpp"
#include "vertex.hpp"

template <int32_t N>
bool slice_cube() {
  const auto vertex_transformations = compute_vertex_transformations(N);
  const auto complexes = compute_cut_complexes_degree_one<N>();
  const auto edges = compute_edges(N);
  const auto edge_transformations =
      compute_edge_transformations(edges, vertex_transformations, N);
  const auto usr = complexes_to_usr<N>(complexes, edges);
  const auto mss = usr_to_mss<N>(usr, edges);
  const auto usr_2 = combine_usr_mss<N>(usr, mss, edges);
  const auto mss_2 = usr_to_mss<N>(usr_2, edges);
  return combine_usr_mss_final<N>(usr, mss);
}

int main() { std::cout << slice_cube<5>() << std::endl; }

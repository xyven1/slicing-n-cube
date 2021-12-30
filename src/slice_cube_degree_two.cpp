#include <cstdint>
#include <iostream>
#include <string>

#include "common.hpp"
#include "edge.hpp"
#include "sliceable_set.hpp"
#include "symmetry.hpp"

int main() {
  constexpr int32_t N = 5;
  const auto vertex_transformations = compute_vertex_transformations(N);
  const auto edges = compute_edges(N);
  const auto edge_transformations =
      compute_edge_transformations(edges, vertex_transformations, N);
  const auto usr = read_from_file<N>(NCUBE_DIR "degree_two/" +
                                     std::to_string(N) + "_usr_1.bin");
  const auto mss = read_from_file<N>(NCUBE_DIR "degree_two/" +
                                     std::to_string(N) + "_usr_1.bin");
  std::cout << usr.size() << std::endl;
  std::cout << mss.size() << std::endl;
}

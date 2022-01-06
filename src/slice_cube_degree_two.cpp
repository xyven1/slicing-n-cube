#include <cstdint>
#include <iostream>
#include <string>

#include "edge.hpp"
#include "sliceable_set.hpp"
#include "vertex.hpp"

using namespace ncube;

int main() {
  constexpr int32_t N = 5;
  const auto edges = compute_edges<N>();
  const auto usr_1 = read_from_file<N>(N_CUBE_OUT_DIR "/degree_two/" +
                                       std::to_string(N) + "_usr_1.bin");
  const auto mss_1 = read_from_file<N>(N_CUBE_OUT_DIR "/degree_two/" +
                                       std::to_string(N) + "_mss_1.bin");
  const auto usr_2 = pairwise_unions<N>(usr_1, mss_1, edges);
  std::cout << "|usr_2| = " << usr_2.size() << std::endl;
  const auto mss_2 = expand_usr<N>(usr_2, edges);
  std::cout << "|mss_2| = " << mss_2.size() << std::endl;
}

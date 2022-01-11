#include <cstdint>
#include <filesystem>
#include <functional>
#include <iostream>
#include <string>

#include "complex.hpp"
#include "edge.hpp"
#include "sliceable_set.hpp"
#include "vertex.hpp"

using namespace ncube;

template <int32_t N>
std::vector<int64_t> edge_cardinalities(
    const std::vector<sliceable_set_t<N>>& sets) {
  std::vector<int64_t> cardinalities(num_edges(N) + 1);
  for (const auto& ss : sets) {
    cardinalities[ss.count()] += 1;
  }
  return cardinalities;
}

template <int32_t N>
std::vector<int64_t> edge_cardinalities_mss(
    std::function<bool(const complex_t<N>&)> is_complex) {
  const auto edges = compute_edges<N>();
  const auto complexes = compute_complexes<N>(is_complex);
  const auto usr = complexes_to_usr<N>(complexes, edges);
  const auto mss = expand_usr<N>(usr, edges);
  return edge_cardinalities<N>(mss);
}

template <int32_t N>
void write_to_file(std::function<bool(const complex_t<N>&)> is_complex,
                   const std::filesystem::path& out) {
  const auto cardinalities = edge_cardinalities_mss<N>(is_complex);
  std::ofstream file(out);
  for (const auto& c : cardinalities) {
    file << c << std::endl;
  }
}

int main() {
  const std::string dir = N_CUBE_OUT_DIR "/cardinality";
  std::filesystem::create_directories(dir);
  write_to_file<2>(is_complex_degree_one<2>, dir + "/degree_one_2.txt");
  write_to_file<3>(is_complex_degree_one<3>, dir + "/degree_one_3.txt");
  write_to_file<4>(is_complex_degree_one<4>, dir + "/degree_one_4.txt");
  write_to_file<5>(is_complex_degree_one<5>, dir + "/degree_one_5.txt");
  write_to_file<6>(is_complex_degree_one<6>, dir + "/degree_one_6.txt");

  write_to_file<2>(is_complex_degree_two<2>, dir + "/degree_two_2.txt");
  write_to_file<3>(is_complex_degree_two<3>, dir + "/degree_two_3.txt");
  write_to_file<4>(is_complex_degree_two<4>, dir + "/degree_two_4.txt");
  write_to_file<5>(is_complex_degree_two<5>, dir + "/degree_two_5.txt");
}

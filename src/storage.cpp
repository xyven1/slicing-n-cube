#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <string>

#include "complex.hpp"
#include "edge.hpp"
#include "sliceable_set.hpp"
#include "vertex.hpp"

using namespace ncube;

template <int32_t N>
void write_degree_two_1_sliceable_sets() {
  const auto edges = compute_edges<N>();
  const auto complexes = compute_complexes<N>(is_complex_degree_two<N>);
  const auto usr = complexes_to_usr<N>(complexes, edges);
  const auto mss = expand_usr<N>(usr, edges);
  std::cout << "degree two |" << N << "_usr_1| = " << usr.size() << std::endl;
  std::cout << "degree two |" << N << "_mss_1| = " << mss.size() << std::endl;
  constexpr auto dir = N_CUBE_OUT_DIR "/degree_two";
  std::filesystem::create_directories(dir);
  const auto path_usr = dir + ("/" + std::to_string(N) + "_usr_1.bin");
  const auto path_mss = dir + ("/" + std::to_string(N) + "_mss_1.bin");
  write_to_file<N>(usr, path_usr);
  write_to_file<N>(mss, path_mss);
}

template <int32_t N>
void write_degree_one_1_sliceable_sets() {
  const auto edges = compute_edges<N>();
  const auto complexes = compute_complexes<N>(is_complex_degree_one<N>);
  const auto usr = complexes_to_usr<N>(complexes, edges);
  const auto mss = expand_usr<N>(usr, edges);
  std::cout << "degree one |" << N << "_usr_1| = " << usr.size() << std::endl;
  std::cout << "degree one |" << N << "_mss_1| = " << mss.size() << std::endl;
  constexpr auto dir = N_CUBE_OUT_DIR "/degree_one";
  std::filesystem::create_directories(dir);
  const auto path_usr = dir + ("/" + std::to_string(N) + "_usr_1.bin");
  const auto path_mss = dir + ("/" + std::to_string(N) + "_mss_1.bin");
  write_to_file<N>(usr, path_usr);
  write_to_file<N>(mss, path_mss);
}

template <int32_t N>
void write_degree_one_1_sliceable_sets_only_usr() {
  const auto edges = compute_edges<N>();
  const auto complexes = compute_complexes<N>(is_complex_degree_one<N>);
  const auto usr = complexes_to_usr<N>(complexes, edges);
  std::cout << "degree one |" << N << "_usr_1| = " << usr.size() << std::endl;
  constexpr auto dir = N_CUBE_OUT_DIR "/degree_one";
  std::filesystem::create_directories(dir);
  const auto path_usr = dir + ("/" + std::to_string(N) + "_usr_1.bin");
  write_to_file<N>(usr, path_usr);
}

template <int32_t N>
void write_degree_one_2_sliceable_sets() {
  const auto edges = compute_edges<N>();
  const auto complexes = compute_complexes<N>(is_complex_degree_one<N>);
  const auto usr_1 = complexes_to_usr<N>(complexes, edges);
  const auto mss_1 = expand_usr<N>(usr_1, edges);
  const auto usr_2 = pairwise_unions<N>(usr_1, mss_1, edges);
  const auto mss_2 = expand_usr<N>(usr_2, edges);
  std::cout << "degree one |" << N << "_usr_1| = " << usr_1.size() << std::endl;
  std::cout << "degree one |" << N << "_mss_1| = " << mss_1.size() << std::endl;
  std::cout << "degree one |" << N << "_usr_2| = " << usr_2.size() << std::endl;
  std::cout << "degree one |" << N << "_mss_2| = " << mss_2.size() << std::endl;
  constexpr auto dir = N_CUBE_OUT_DIR "/degree_one";
  std::filesystem::create_directories(dir);
  const auto path_usr_1 = dir + ("/" + std::to_string(N) + "_usr_1.bin");
  const auto path_mss_1 = dir + ("/" + std::to_string(N) + "_mss_1.bin");
  const auto path_usr_2 = dir + ("/" + std::to_string(N) + "_usr_2.bin");
  const auto path_mss_2 = dir + ("/" + std::to_string(N) + "_mss_2.bin");
  write_to_file<N>(usr_1, path_usr_1);
  write_to_file<N>(mss_1, path_mss_1);
  write_to_file<N>(usr_2, path_usr_2);
  write_to_file<N>(mss_2, path_mss_2);
}

int main() {
  write_degree_one_2_sliceable_sets<2>();
  write_degree_one_2_sliceable_sets<3>();
  write_degree_one_2_sliceable_sets<4>();
  write_degree_one_2_sliceable_sets<5>();
  write_degree_one_1_sliceable_sets<6>();
  write_degree_one_1_sliceable_sets_only_usr<7>();

  write_degree_two_1_sliceable_sets<2>();
  write_degree_two_1_sliceable_sets<3>();
  write_degree_two_1_sliceable_sets<4>();
  write_degree_two_1_sliceable_sets<5>();
}

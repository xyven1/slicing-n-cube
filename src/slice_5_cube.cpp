#include <chrono>
#include <iostream>

#include "sliceable_set.hpp"

using namespace ncube;

int main() {
  constexpr auto usr_2_path = N_CUBE_OUT_DIR "/degree_one/5_usr_2.bin";
  constexpr auto mss_2_path = N_CUBE_OUT_DIR "/degree_one/5_mss_2.bin";
  const auto usr_2 = read_from_file<5>(usr_2_path);
  const auto mss_2 = read_from_file<5>(mss_2_path);
  const auto start = std::chrono::high_resolution_clock::now();
  const auto slices_cube = pairwise_unions_slice_cube<5>(usr_2, mss_2);
  const auto stop = std::chrono::high_resolution_clock::now();
  const auto duration =
      std::chrono::duration_cast<std::chrono::seconds>(stop - start);
  std::cout << "Execution time of pairwise_unions_slice_cube: "
            << duration.count() << std::endl;
  std::cout << "Can four hyperplanes slice the 5-cube: " << slices_cube
            << std::endl;
}

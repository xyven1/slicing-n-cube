#include <cstdint>
#include <iostream>

#include "edge.hpp"
#include "low_weight.hpp"
#include "vertex.hpp"

using namespace ncube;

template <int32_t N>
int32_t slice_cube_one_weight(const std::vector<int32_t>& thresholds) {
  const auto edges = compute_edges<N>();
  const auto sets = compute_one_weight_mss<N>(thresholds, edges);
  const auto k = slice_cube_low_weight<N>(sets, edges);
  return k;
}

int main() {
  std::vector<int32_t> thresholds = {0, 1};
  std::cout << "Minimum number of halfspaces with normal vector in {-1, 1} and "
               "threshold in {0, 1} required to slice n-cube"
            << std::endl;
  std::cout << "n = 2: " << slice_cube_one_weight<2>(thresholds) << std::endl;
  std::cout << "n = 3: " << slice_cube_one_weight<3>(thresholds) << std::endl;
  std::cout << "n = 4: " << slice_cube_one_weight<4>(thresholds) << std::endl;
  std::cout << "n = 5: " << slice_cube_one_weight<5>(thresholds) << std::endl;
  std::cout << "n = 6: " << slice_cube_one_weight<6>(thresholds) << std::endl;
  std::cout << "n = 7: " << slice_cube_one_weight<7>(thresholds) << std::endl;
  std::cout << "Minimum number of halfspaces with normal vector in {-1, 1} and "
               "threshold in {0, ..., n} required to slice n-cube"
            << std::endl;
  thresholds.push_back(2);
  std::cout << "n = 2: " << slice_cube_one_weight<2>(thresholds) << std::endl;
  thresholds.push_back(3);
  std::cout << "n = 3: " << slice_cube_one_weight<3>(thresholds) << std::endl;
  thresholds.push_back(4);
  std::cout << "n = 4: " << slice_cube_one_weight<4>(thresholds) << std::endl;
  thresholds.push_back(5);
  std::cout << "n = 5: " << slice_cube_one_weight<5>(thresholds) << std::endl;
  thresholds.push_back(6);
  std::cout << "n = 6: " << slice_cube_one_weight<6>(thresholds) << std::endl;
}

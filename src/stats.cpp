#include <cstdint>
#include <iostream>

#include "complex.hpp"
#include "edge.hpp"
#include "sliceable_set.hpp"
#include "vertex.hpp"

using namespace ncube;

template <int32_t N>
std::vector<int64_t> compute_edge_frequencies() {
  const auto edges = compute_edges<N>();
  const auto complexes = compute_complexes<N>(is_complex_degree_one<N>);
  const auto usr = complexes_to_usr<N>(complexes, edges);
  const auto mss = expand_usr<N>(usr, edges);
  std::vector<int64_t> frequencies(num_edges(N) + 1);
  for (const auto& ss : mss) {
    frequencies[ss.count()] += 1;
  }
  return frequencies;
}

template <int32_t N>
int64_t count_pairs_cardinality(const std::vector<sliceable_set_t<N>>& usr,
                                const std::vector<sliceable_set_t<N>>& mss) {
  constexpr auto max_cardinality = sliceable_set_t<N>().size();
  std::vector<int64_t> cardinalities(max_cardinality + 1);
  for (const auto& ss : mss) {
    cardinalities[ss.count()] += 1;
  }
  int64_t count = 0;
  for (const auto& ss : usr) {
    for (auto i = max_cardinality - ss.count(); i < cardinalities.size(); ++i) {
      count += cardinalities[i];
    }
  }
  return count;
}

template <int32_t N>
int64_t count_pairs_leading_zeros(const std::vector<sliceable_set_t<N>>& usr,
                                  const std::vector<sliceable_set_t<N>>& mss) {
  std::vector<int64_t> leading_ones(sliceable_set_t<N>().size() + 1);
  for (const auto& ss : mss) {
    leading_ones[get_leading_ones<N>(ss)] += 1;
  }
  int64_t count = 0;
  for (const auto& ss : usr) {
    for (auto i = get_leading_zeros<N>(ss); i < leading_ones.size(); ++i) {
      count += leading_ones[i];
    }
  }
  return count;
}

int main() {
  for (const auto& freq : compute_edge_frequencies<5>()) {
    std::cout << freq << std::endl;
  }
}

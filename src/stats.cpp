#include <cstdint>
#include <iostream>

#include "common.hpp"
#include "complex.hpp"
#include "edge.hpp"
#include "sliceable_set.hpp"
#include "symmetry.hpp"

template <int32_t N>
std::vector<int64_t> compute_edge_frequencies() {
  const std::vector<vertex_trans_t> vertex_transformations =
      compute_vertex_transformations(N);
  const std::vector<complex_t<N>> complexes =
      compute_cut_complexes_degree_one<N>(vertex_transformations);
  const std::vector<edge_t> edges = compute_edges(N);
  const std::vector<sliceable_set_t<N>> mss =
      complexes_to_mss<N>(complexes, vertex_transformations, edges);
  std::vector<int64_t> frequencies(num_edges(N) + 1);
  for (const sliceable_set_t<N>& ss : mss) {
    frequencies[ss.count()] += 1;
  }
  return frequencies;
}

template <int32_t N>
int64_t count_pairs_cardinality(const std::vector<sliceable_set_t<N>>& usr,
                                const std::vector<sliceable_set_t<N>>& mss) {
  constexpr std::size_t max_cardinality = sliceable_set_t<N>().size();
  std::vector<int64_t> cardinalities(max_cardinality + 1);
  for (const sliceable_set_t<N>& ss : mss) {
    cardinalities[ss.count()] += 1;
  }
  int64_t count = 0;
  for (const sliceable_set_t<N>& ss : usr) {
    for (std::size_t i = max_cardinality - ss.count(); i < cardinalities.size();
         ++i) {
      count += cardinalities[i];
    }
  }
  return count;
}

template <int32_t N>
int64_t count_pairs_leading_zeros(const std::vector<sliceable_set_t<N>>& usr,
                                  const std::vector<sliceable_set_t<N>>& mss) {
  std::vector<int64_t> leading_ones(sliceable_set_t<N>().size() + 1);
  for (const sliceable_set_t<N>& ss : mss) {
    leading_ones[get_leading_ones<N>(ss)] += 1;
  }
  int64_t count = 0;
  for (const sliceable_set_t<N>& ss : usr) {
    for (std::size_t i = get_leading_zeros<N>(ss); i < leading_ones.size();
         ++i) {
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

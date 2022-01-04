#include <cstdint>
#include <iostream>

#include "common.hpp"
#include "complex.hpp"
#include "edge.hpp"
#include "sliceable_set.hpp"
#include "symmetry.hpp"

template <int32_t N>
std::vector<int64_t> compute_edge_frequencies() {
  const auto vertex_transformations = compute_vertex_transformations(N);
  const auto complexes =
      compute_cut_complexes_degree_one<N>(vertex_transformations);
  const auto edges = compute_edges(N);
  const auto edge_transformations =
      compute_edge_transformations(edges, vertex_transformations, N);
  const auto usr = complexes_to_usr<N>(complexes, edge_transformations, edges);
  const auto mss = usr_to_mss<N>(usr, edge_transformations);
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
  for (const auto& ss : mss) {
    leading_ones[get_leading_ones<N>(ss)] += 1;
  }
  int64_t count = 0;
  for (const auto& ss : usr) {
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

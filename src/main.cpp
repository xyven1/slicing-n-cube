#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_set>

#include "common.hpp"
#include "complex.hpp"
#include "edge.hpp"
#include "low_weight.hpp"
#include "sliceable_set.hpp"
#include "symmetry.hpp"

template <int32_t N>
void write_two_sliceable_sets() {
  const std::vector<vertex_trans_t> vertex_transformations =
      compute_vertex_transformations(N);
  const std::vector<complex_t<N>> complexes =
      compute_cut_complexes<N>(vertex_transformations);
  const std::vector<edge_t> edges = compute_edges(N);
  const std::vector<edge_trans_t> edge_transformations =
      compute_edge_transformations(edges, vertex_transformations, N);
  const std::vector<sliceable_set_t<N>> usr =
      complexes_to_usr<N>(complexes, edge_transformations, edges);
  const std::vector<sliceable_set_t<N>> mss =
      complexes_to_mss<N>(complexes, vertex_transformations, edges);
  std::cout << usr.size() << std::endl;
  std::cout << mss.size() << std::endl;
  const std::vector<sliceable_set_t<N>> usr_2 =
      combine_usr_mss<N>(usr, mss, edge_transformations);
  std::cout << usr_2.size() << std::endl;
  const std::vector<sliceable_set_t<N>> mss_2 =
      usr_to_mss<N>(usr_2, edge_transformations);
  std::cout << mss_2.size() << std::endl;
  write_to_file<N>(usr, NCUBE_DIR + std::to_string(N) + "_usr_1.bin");
  write_to_file<N>(mss, NCUBE_DIR + std::to_string(N) + "_mss_1.bin");
  write_to_file<N>(usr_2, NCUBE_DIR + std::to_string(N) + "_usr_2.bin");
  write_to_file<N>(mss_2, NCUBE_DIR + std::to_string(N) + "_mss_2.bin");
}

template <int32_t N>
void print_sizes() {
  const std::string usr_path_1 = NCUBE_DIR + std::to_string(N) + "_usr_1.bin";
  const std::string mss_path_1 = NCUBE_DIR + std::to_string(N) + "_mss_1.bin";
  const std::string usr_path_2 = NCUBE_DIR + std::to_string(N) + "_usr_2.bin";
  const std::string mss_path_2 = NCUBE_DIR + std::to_string(N) + "_mss_2.bin";
  const std::vector<sliceable_set_t<N>> usr_1 = read_from_file<N>(usr_path_1);
  const std::vector<sliceable_set_t<N>> mss_1 = read_from_file<N>(mss_path_1);
  const std::vector<sliceable_set_t<N>> usr_2 = read_from_file<N>(usr_path_2);
  const std::vector<sliceable_set_t<N>> mss_2 = read_from_file<N>(mss_path_2);
  std::cout << "n = " << N << std::endl;
  std::cout << "usr_1 size = " << usr_1.size() << std::endl;
  std::cout << "mss_1 size = " << mss_1.size() << std::endl;
  std::cout << "usr_2 size = " << usr_2.size() << std::endl;
  std::cout << "mss_2 size = " << mss_2.size() << std::endl;
}

template <int32_t N>
std::vector<int64_t> compute_edge_frequencies() {
  const std::vector<vertex_trans_t> vertex_transformations =
      compute_vertex_transformations(N);
  const std::vector<complex_t<N>> complexes =
      compute_cut_complexes<N>(vertex_transformations);
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
  constexpr int32_t N = 5;
  const auto edges = compute_edges(N);
  const auto vertex_transformations = compute_vertex_transformations(N);
  const auto edge_transformations =
      compute_edge_transformations(edges, vertex_transformations, N);
  std::vector<double> distances;
  for (double i = -N; i <= N; i += 0.5) {
    distances.push_back(i);
  }
  const auto low_weight_sets =
      compute_low_weight_sliceable_sets<N>(distances, edges);
  const auto k = combine_low_weight_sliceable_sets<N>(low_weight_sets,
                                                      edge_transformations);
  std::cout << k << std::endl;
}

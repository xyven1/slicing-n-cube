#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_set>

#include "common.hpp"
#include "complex.hpp"
#include "edge.hpp"
#include "sliceable_set.hpp"
#include "symmetry.hpp"

void write_two_sliceable_sets(int32_t n) {
  const std::vector<vertex_trans_t> vertex_transformations =
      compute_vertex_transformations(n);
  const std::vector<complex_t> complexes =
      compute_cut_complexes(vertex_transformations, n);
  const std::vector<edge_t> edges = compute_edges(n);
  const std::vector<edge_trans_t> edge_transformations =
      compute_edge_transformations(edges, vertex_transformations, n);
  const std::vector<sliceable_set_t> usr =
      complexes_to_usr(complexes, edge_transformations, edges, n);
  const std::vector<sliceable_set_t> mss =
      complexes_to_mss(complexes, vertex_transformations, edges, n);
  std::cout << usr.size() << std::endl;
  std::cout << mss.size() << std::endl;
  const std::vector<sliceable_set_t> usr_2 =
      combine_usr_mss(usr, mss, edge_transformations, n);
  std::cout << usr_2.size() << std::endl;
  const std::vector<sliceable_set_t> mss_2 =
      usr_to_mss(usr_2, edge_transformations, n);
  std::cout << mss_2.size() << std::endl;
  write_to_file(usr, NCUBE_DIR + std::to_string(n) + "_usr_1.bin");
  write_to_file(mss, NCUBE_DIR + std::to_string(n) + "_mss_1.bin");
  write_to_file(usr_2, NCUBE_DIR + std::to_string(n) + "_usr_2.bin");
  write_to_file(mss_2, NCUBE_DIR + std::to_string(n) + "_mss_2.bin");
}

void print_sizes(int32_t n) {
  const std::string usr_path_1 = NCUBE_DIR + std::to_string(n) + "_usr_1.bin";
  const std::string mss_path_1 = NCUBE_DIR + std::to_string(n) + "_mss_1.bin";
  const std::string usr_path_2 = NCUBE_DIR + std::to_string(n) + "_usr_2.bin";
  const std::string mss_path_2 = NCUBE_DIR + std::to_string(n) + "_mss_2.bin";
  const std::vector<sliceable_set_t> usr_1 = read_from_file(usr_path_1);
  const std::vector<sliceable_set_t> mss_1 = read_from_file(mss_path_1);
  const std::vector<sliceable_set_t> usr_2 = read_from_file(usr_path_2);
  const std::vector<sliceable_set_t> mss_2 = read_from_file(mss_path_2);
  std::cout << "n = " << n << std::endl;
  std::cout << "usr_1 size = " << usr_1.size() << std::endl;
  std::cout << "mss_1 size = " << mss_1.size() << std::endl;
  std::cout << "usr_2 size = " << usr_2.size() << std::endl;
  std::cout << "mss_2 size = " << mss_2.size() << std::endl;
}

std::vector<int64_t> compute_edge_frequencies(int32_t n) {
  const std::vector<vertex_trans_t> vertex_transformations =
      compute_vertex_transformations(n);
  const std::vector<complex_t> complexes =
      compute_cut_complexes(vertex_transformations, n);
  const std::vector<edge_t> edges = compute_edges(n);
  const std::vector<sliceable_set_t> mss =
      complexes_to_mss(complexes, vertex_transformations, edges, n);
  std::vector<int64_t> frequencies(num_edges(n) + 1);
  for (const sliceable_set_t& ss : mss) {
    frequencies[ss.count()] += 1;
  }
  return frequencies;
}

int64_t count_pairs_cardinality(const std::vector<sliceable_set_t>& usr,
                                const std::vector<sliceable_set_t>& mss) {
  constexpr std::size_t max_cardinality = sliceable_set_t().size();
  std::vector<int64_t> cardinalities(max_cardinality + 1);
  for (const sliceable_set_t& ss : mss) {
    cardinalities[ss.count()] += 1;
  }
  int64_t count = 0;
  for (const sliceable_set_t& ss : usr) {
    for (std::size_t i = max_cardinality - ss.count(); i < cardinalities.size();
         ++i) {
      count += cardinalities[i];
    }
  }
  return count;
}

int64_t count_pairs_leading_zeros(const std::vector<sliceable_set_t>& usr,
                                  const std::vector<sliceable_set_t>& mss) {
  std::vector<int64_t> leading_ones(sliceable_set_t().size() + 1);
  for (const sliceable_set_t& ss : mss) {
    leading_ones[get_leading_ones(ss)] += 1;
  }
  int64_t count = 0;
  for (const sliceable_set_t& ss : usr) {
    for (std::size_t i = get_leading_zeros(ss); i < leading_ones.size(); ++i) {
      count += leading_ones[i];
    }
  }
  return count;
}

int main() {
  constexpr int32_t n = N;
  const std::string usr_path = NCUBE_DIR + std::to_string(n) + "_usr_2.bin";
  const std::string mss_path = NCUBE_DIR + std::to_string(n) + "_mss_2.bin";
  const std::vector<sliceable_set_t> usr = read_from_file(usr_path);
  const std::vector<sliceable_set_t> mss = read_from_file(mss_path);
  const bool slices_all = combine_usr_mss_final(usr, mss);
  std::cout << slices_all << std::endl;
}

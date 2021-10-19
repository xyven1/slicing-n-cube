#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_set>

#include "common.hpp"
#include "complex.hpp"
#include "edge.hpp"
#include "sliceable_set.hpp"
#include "symmetry.hpp"

std::string vertex_to_str(vertex_t v, int32_t n) {
  std::string str;
  str.reserve(n);
  for (int32_t i = 0; i < n; ++i) {
    const int32_t bit_i = (v >> i) & 1;
    str += (bit_i) ? "1" : "0";
  }
  return str;
}

void evaluate_f() {
  const std::vector<double> f = {1, 1, 1, 0.9};
  const int32_t n = static_cast<int32_t>(f.size() - 1);
  for (vertex_t v = 0; v < num_vertices(n); ++v) {
    double x = f[n];
    for (int32_t i = 0; i < n; ++i) {
      const int32_t val_i = (v >> i) & 1;
      const double coordinate_i = (val_i) ? 1 : -1;
      x += f[i] * coordinate_i;
    }
    std::cout << "vertex " << v << ": " << x << std::endl;
  }
}

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

int main() {
  constexpr int32_t n = N;
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
  const std::vector<sliceable_set_t> usr_2 =
      combine_usr_mss(usr, mss, edge_transformations, n);
  const std::vector<sliceable_set_t> mss_2 =
      usr_to_mss(usr_2, edge_transformations, n);
  if (usr != read_from_file(NCUBE_DIR + std::to_string(n) + "_usr_1.bin")) {
    throw std::runtime_error("usr");
  }
  if (mss != read_from_file(NCUBE_DIR + std::to_string(n) + "_mss_1.bin")) {
    throw std::runtime_error("mss");
  }
  if (usr_2 != read_from_file(NCUBE_DIR + std::to_string(n) + "_usr_2.bin")) {
    throw std::runtime_error("usr_2");
  }
  if (mss_2 != read_from_file(NCUBE_DIR + std::to_string(n) + "_mss_2.bin")) {
    throw std::runtime_error("mss_2");
  }
}

#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_set>

#include "common.hpp"
#include "complex.hpp"
#include "edge.hpp"
#include "symmetry.hpp"

std::string vertex_to_str(vertex_t v, int32_t n) {
  std::string str;
  str.reserve(n);
  for (int32_t i = 0; i < n; ++i) {
    const uint32_t bit_i = (v >> i) & 1;
    str += (bit_i) ? "1" : "0";
  }
  return str;
}

void evaluate_f() {
  const std::vector<double> f = {1, 1, 1, 0.9};
  const int32_t n = static_cast<int32_t>(f.size() - 1);
  for (vertex_t v = 0; v < (1u << n); ++v) {
    double x = f[n];
    for (int32_t i = 0; i < n; ++i) {
      const uint32_t val_i = (v >> i) & 1;
      const double coordinate_i = (val_i) ? 1 : -1;
      x += f[i] * coordinate_i;
    }
    std::cout << "vertex " << v << ": " << x << std::endl;
  }
}

int main() {
  constexpr int32_t n = 6;
  static_assert((1 << n) == complex_t().size());
  static_assert(n * (1 << (n - 1)) == sliceable_set_t().size());
  const std::vector<symmetry_t> symmetries = compute_symmetries(n);
  const std::vector<inversion_t> inversions =
      compute_vertex_inversions(symmetries, n);
  const std::vector<transformation_t> transformations =
      compute_vertex_transformations(symmetries, n);
  const std::vector<complex_t> complexes = compute_cut_complexes(inversions, n);
  std::cout << complexes.size() << std::endl;
  std::unordered_set<sliceable_set_t> mss;
  const std::vector<edge_t> edges = compute_edges(n);
  for (const complex_t& complex : complexes) {
    for (const complex_t& expansion : expand_complex(transformations, complex, n)) {
      mss.insert(complex_to_sliceable_set(expansion, edges, n));
    }
  }
  std::cout << mss.size() << std::endl;
}

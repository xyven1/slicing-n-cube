#include <cstdint>
#include <iostream>
#include <string>

#include "common.hpp"
#include "complex.hpp"

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
  const std::vector<complex_t> complexes = compute_cut_complexes(n);
  std::cout << complexes.size() << std::endl;
  // for (const complex_t& c : complexes) {
  //   std::cout << complex_to_edges(c, n) << std::endl;
  // }
}
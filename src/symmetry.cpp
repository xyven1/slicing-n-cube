#include "symmetry.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

#include "common.hpp"

std::vector<symmetry_t> compute_symmetries(int32_t n) {
  const int32_t two_to_n = pow(2, n);
  const int32_t size = two_to_n * factorial(n);
  std::vector<symmetry_t> symmetries;
  symmetries.reserve(size);

  std::vector<int32_t> permutation(n);
  for (int32_t i = 0; i < n; ++i) {
    permutation[i] = i;
  }
  do {
    // if the i-th least significant bit is set, negate the i-th coordinate
    for (uint32_t signs = 0; signs < (1u << n); ++signs) {
      symmetry_t symmetry;
      symmetry.reserve(n);
      for (int32_t i = 0; i < n; ++i) {
        const uint32_t mask = 1 << i;
        const int32_t sign = (mask & signs) ? 1 : 0;
        symmetry.emplace_back(sign, permutation[i]);
      }
      symmetries.push_back(symmetry);
    }
  } while (std::next_permutation(permutation.begin(), permutation.end()));
  return symmetries;
}

vertex_t transform_vertex(const symmetry_t& sym, vertex_t v, int32_t n) {
  vertex_t transformation = 0;
  for (int32_t i = 0; i < n; ++i) {
    const uint32_t val_i = (v >> i) & 1;
    const uint32_t val_i_sign = val_i ^ sym[i].first;
    const uint32_t val_i_sign_and_position = val_i_sign << sym[i].second;
    transformation |= val_i_sign_and_position;
  }
  return transformation;
}

vertex_t transform_vertex_inv(const symmetry_t& sym, vertex_t v, int32_t n) {
  vertex_t transformation = 0;
  for (int32_t i = 0; i < n; ++i) {
    const uint32_t val_i = (v >> sym[i].second) & 1;
    const uint32_t val_i_sign = val_i ^ sym[i].first;
    const uint32_t val_i_sign_and_position = val_i_sign << i;
    transformation |= val_i_sign_and_position;
  }
  return transformation;
}

complex_t transform_complex(const complex_t& complex, const symmetry_t& sym,
                            int32_t n) {
  // std::cout << "Transforming complex " << complex << " with " << sym <<
  // std::endl;
  complex_t transformation;
  for (vertex_t v = 0; v < (1u << n); ++v) {
    if (complex[v]) {
      transformation[transform_vertex(sym, v, n)] = true;
    }
  }
  // std::cout << "Transformation is " << transformation << std::endl;
  return transformation;
}

complex_t transform_complex_and_min(const complex_t& complex,
                                    const symmetry_t& sym, int32_t n,
                                    const complex_t& min_complex) {
  complex_t transformation;
  for (vertex_t v = (1u << n) - 1; v < (1u << n); --v) {
    const vertex_t inv = transform_vertex_inv(sym, v, n);
    transformation[v] = complex[inv];
    if (transformation[v] && !min_complex[v]) {
      return min_complex;
    }
    if (!transformation[v] && min_complex[v]) {
      return transform_complex(complex, sym, n);
    }
  }
  return transformation;
}

#include <cstdint>
#include <vector>

#include "common.hpp"
#include "complex.hpp"
#include "lp.hpp"
#include "symmetry.hpp"

complex_t unique_complex(const complex_t& complex,
                         const std::vector<symmetry_t>& symmetries, int32_t n) {
  complex_t min_transformation(complex);
  for (const symmetry_t& sym : symmetries) {
    const complex_t transformation = transform_complex(complex, sym, n);
    if (min_transformation > transformation) {
      min_transformation = transformation;
    }
  }
  return min_transformation;
}

std::vector<vertex_t> adjacent_vertices_of_complex(const complex_t& complex,
                                                   int32_t n) {
  // TODO: Use set
  std::vector<vertex_t> adjacent_vertices;
  for (const vertex_t& v : complex) {
    for (int32_t i = 0; i < n; ++i) {
      vertex_t neighbour = v ^ (1 << i);
      if (complex.count(neighbour) == 0 &&
          std::find(adjacent_vertices.begin(), adjacent_vertices.end(),
                    neighbour) == adjacent_vertices.end()) {
        adjacent_vertices.push_back(neighbour);
      }
    }
  }
  return adjacent_vertices;
}

std::vector<edge_t> complex_to_edges(const complex_t& complex, int32_t n) {
  std::vector<edge_t> edges;
  for (const vertex_t& v : complex) {
    for (int32_t i = 0; i < n; ++i) {
      const vertex_t neighbour = v ^ (1 << i);
      if (complex.count(neighbour) == 0) {
        if (v < neighbour) {
          edges.emplace_back(v, neighbour);
        } else {
          edges.emplace_back(neighbour, v);
        }
      }
    }
  }
  return edges;
}

std::vector<complex_t> compute_cut_complexes(int32_t n) {
  std::vector<symmetry_t> symmetries = compute_symmetries(n);
  const int32_t l = pow(2, n - 1);
  // There is exactly one USR of a cut complex of size l=1
  std::vector<complex_t> complexes = {{0}};
  std::size_t prev_begin = 0;
  std::size_t prev_end = complexes.size();
  for (int32_t i = 1; i < l; ++i) {
    // std::cout << "Expanding complexes of size " << i << std::endl;
    for (std::size_t j = prev_begin; j < prev_end; ++j) {
      // std::cout << "Expanding complex " << complexes[j] << std::endl;
      for (const vertex_t& v : adjacent_vertices_of_complex(complexes[j], n)) {
        // std::cout << "\nTrying adjacent vertex " << v << std::endl;
        complex_t maybe_complex(complexes[j]);
        maybe_complex.insert(v);
        maybe_complex = unique_complex(maybe_complex, symmetries, n);
        // std::cout << "Maybe new complex: " << maybe_complex << std::endl;
        if (std::find(complexes.begin() + prev_end, complexes.end(),
                      maybe_complex) == complexes.end()) {
          if (is_complex(maybe_complex, n)) {
            complexes.push_back(maybe_complex);
            // std::cout << "New complex: " << maybe_complex << std::endl;
          }
        } else {
          // std::cout << "Already known complex" << std::endl;
        }
      }
    }
    prev_begin = prev_end;
    prev_end = complexes.size();
  }
  return complexes;
}

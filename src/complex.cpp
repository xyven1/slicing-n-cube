#include "complex.hpp"

#include <cstdint>
#include <vector>

#include "common.hpp"
#include "edge.hpp"
#include "lp.hpp"
#include "symmetry.hpp"

complex_t unique_complex(const complex_t& complex,
                         const std::vector<std::vector<vertex_t>>& inversions,
                         int32_t n) {
  complex_t min_complex(complex);
  for (const std::vector<vertex_t>& inversion : inversions) {
    complex_t transformation;
    bool is_new_min = false;
    for (vertex_t v = (1u << n) - 1; v < (1u << n); --v) {
      const vertex_t inv = inversion[v];
      transformation[v] = complex[inv];
      if (!is_new_min && transformation[v] && !min_complex[v]) {
        // min_complex is still smaller
        break;
      }
      is_new_min |= !transformation[v] && min_complex[v];
    }
    if (is_new_min) {
      min_complex = transformation;
    }
  }
  return min_complex;
}

std::vector<vertex_t> adjacent_vertices_of_complex(const complex_t& complex,
                                                   int32_t n) {
  std::vector<vertex_t> adjacent_vertices;
  for (vertex_t v = 0; v < (1u << n); ++v) {
    if (complex[v]) {
      for (int32_t i = 0; i < n; ++i) {
        vertex_t neighbour = v ^ (1 << i);
        if (!complex[neighbour] &&
            std::find(adjacent_vertices.begin(), adjacent_vertices.end(),
                      neighbour) == adjacent_vertices.end()) {
          adjacent_vertices.push_back(neighbour);
        }
      }
    }
  }
  return adjacent_vertices;
}

sliceable_set_t complex_to_sliceable_set(const complex_t& complex,
                                         const std::vector<edge_t>& edges,
                                         int32_t n) {
  sliceable_set_t sliceable_set;
  for (vertex_t v = 0; v < (1u << n); ++v) {
    if (complex[v]) {
      for (int32_t i = 0; i < n; ++i) {
        const vertex_t neighbour = v ^ (1 << i);
        if (!complex[neighbour]) {
          const edge_t e =
              (v < neighbour) ? edge_t(v, neighbour) : edge_t(neighbour, v);
          sliceable_set[edge_to_int(e, edges)] = true;
        }
      }
    }
  }
  return sliceable_set;
}

std::vector<complex_t> compute_cut_complexes(int32_t n) {
  const std::vector<symmetry_t> symmetries = compute_symmetries(n);
  const std::vector<std::vector<vertex_t>> inversions =
      compute_vertex_inversions(symmetries, n);
  const int32_t l = pow(2, n - 1);
  // There is exactly one USR of a cut complex of size l=1
  std::vector<complex_t> complexes = {{1}};
  std::size_t prev_begin = 0;
  std::size_t prev_end = complexes.size();
  for (int32_t i = 1; i < l; ++i) {
    // std::cout << "Expanding complexes of size " << i << std::endl;
    for (std::size_t j = prev_begin; j < prev_end; ++j) {
      // std::cout << "Expanding complex " << complexes[j] << std::endl;
      for (const vertex_t& v : adjacent_vertices_of_complex(complexes[j], n)) {
        // std::cout << "\nTrying adjacent vertex " << v << std::endl;
        complex_t maybe_complex(complexes[j]);
        maybe_complex[v] = true;
        maybe_complex = unique_complex(maybe_complex, inversions, n);
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

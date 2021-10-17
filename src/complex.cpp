#include "complex.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

#include "common.hpp"
#include "edge.hpp"
#include "lp.hpp"
#include "sliceable_set.hpp"
#include "symmetry.hpp"

complex_t unique_complex(const complex_t& complex,
                         const std::vector<vertex_trans_t>& transformations,
                         int32_t n) {
  complex_t min_complex(complex);
  for (const vertex_trans_t& transformation : transformations) {
    complex_t complex_trans;
    bool is_new_min = false;
    for (vertex_t v = num_vertices(n) - 1; v >= 0; --v) {
      const vertex_t v_inversion = transformation[v];
      complex_trans[v] = complex[v_inversion];
      if (!is_new_min && complex_trans[v] && !min_complex[v]) {
        // min_complex is still smaller
        break;
      }
      is_new_min |= !complex_trans[v] && min_complex[v];
    }
    if (is_new_min) {
      min_complex = complex_trans;
    }
  }
  return min_complex;
}

std::vector<vertex_t> adjacent_vertices_of_complex(const complex_t& complex,
                                                   int32_t n) {
  std::vector<vertex_t> adjacent_vertices;
  for (vertex_t v = 0; v < num_vertices(n); ++v) {
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
  for (vertex_t v = 0; v < num_vertices(n); ++v) {
    if (complex[v]) {
      for (int32_t i = 0; i < n; ++i) {
        const vertex_t u = v ^ (1 << i);
        if (!complex[u]) {
          const edge_t e = (u < v) ? edge_t(u, v) : edge_t(v, u);
          sliceable_set[edge_to_int(e, edges)] = true;
        }
      }
    }
  }
  return sliceable_set;
}

std::vector<sliceable_set_t> complexes_to_usr(
    const std::vector<complex_t>& complexes,
    const std::vector<edge_trans_t>& edge_transformations,
    const std::vector<edge_t>& edges, int32_t n) {
  std::vector<sliceable_set_t> usr;
  for (const complex_t& complex : complexes) {
    sliceable_set_t ss = complex_to_sliceable_set(complex, edges, n);
    ss = unique_sliceable_set(ss, edge_transformations, n);
    usr.push_back(ss);
  }
  return usr;
}

std::vector<sliceable_set_t> complexes_to_mss(
    const std::vector<complex_t>& complexes,
    const std::vector<vertex_trans_t>& vertex_transformations,
    const std::vector<edge_t>& edges, int32_t n) {
  std::unordered_set<sliceable_set_t> mss;
  for (const complex_t& complex : complexes) {
    for (const vertex_trans_t& transformation : vertex_transformations) {
      complex_t complex_trans;
      for (vertex_t v = 0; v < num_vertices(n); ++v) {
        complex_trans[v] = complex[transformation[v]];
      }
      mss.insert(complex_to_sliceable_set(complex_trans, edges, n));
    }
  }
  return std::vector<sliceable_set_t>(mss.begin(), mss.end());
}

std::vector<complex_t> compute_cut_complexes(
    const std::vector<vertex_trans_t>& transformations, int32_t n) {
  const int32_t l = num_vertices(n) / 2;
  // There is exactly one USR of all complexes of size 1
  std::vector<complex_t> complexes = {{1}};
  // The range [prev_begin, prev_end) contains all complexes of size i
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
        maybe_complex = unique_complex(maybe_complex, transformations, n);
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

#include "complex.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

#include "common.hpp"
#include "edge.hpp"
#include "lp.hpp"
#include "sliceable_set.hpp"
#include "symmetry.hpp"

template <int32_t N>
complex_t<N> unique_complex(
    const complex_t<N>& complex,
    const std::vector<vertex_trans_t>& transformations) {
  complex_t<N> min_complex(complex);
  for (const vertex_trans_t& transformation : transformations) {
    complex_t<N> complex_trans;
    bool is_new_min = false;
    for (vertex_t v = num_vertices(N) - 1; v >= 0; --v) {
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

template <int32_t N>
std::vector<vertex_t> adjacent_vertices_of_complex(
    const complex_t<N>& complex) {
  std::vector<vertex_t> adjacent_vertices;
  for (vertex_t v = 0; v < num_vertices(N); ++v) {
    if (complex[v]) {
      for (int32_t i = 0; i < N; ++i) {
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

template <int32_t N>
sliceable_set_t<N> complex_to_sliceable_set(const complex_t<N>& complex,
                                            const std::vector<edge_t>& edges) {
  sliceable_set_t<N> sliceable_set;
  for (vertex_t v = 0; v < num_vertices(N); ++v) {
    if (complex[v]) {
      for (int32_t i = 0; i < N; ++i) {
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

template <int32_t N>
std::vector<sliceable_set_t<N>> complexes_to_usr(
    const std::vector<complex_t<N>>& complexes,
    const std::vector<edge_trans_t>& edge_transformations,
    const std::vector<edge_t>& edges) {
  std::vector<sliceable_set_t<N>> usr;
  for (const complex_t<N>& complex : complexes) {
    sliceable_set_t<N> ss = complex_to_sliceable_set<N>(complex, edges);
    ss = unique_sliceable_set<N>(ss, edge_transformations);
    usr.push_back(ss);
  }
  return usr;
}

template <int32_t N>
std::vector<sliceable_set_t<N>> complexes_to_mss(
    const std::vector<complex_t<N>>& complexes,
    const std::vector<vertex_trans_t>& vertex_transformations,
    const std::vector<edge_t>& edges) {
  std::unordered_set<sliceable_set_t<N>> mss;
  for (const complex_t<N>& complex : complexes) {
    for (const vertex_trans_t& transformation : vertex_transformations) {
      complex_t<N> complex_trans;
      for (vertex_t v = 0; v < num_vertices(N); ++v) {
        complex_trans[v] = complex[transformation[v]];
      }
      mss.insert(complex_to_sliceable_set<N>(complex_trans, edges));
    }
  }
  return std::vector<sliceable_set_t<N>>(mss.begin(), mss.end());
}

template <int32_t N>
std::vector<complex_t<N>> compute_cut_complexes(
    const std::vector<vertex_trans_t>& transformations) {
  const int32_t l = num_vertices(N) / 2;
  // There is exactly one USR of all complexes of size 1
  std::vector<complex_t<N>> complexes = {{1}};
  // The range [prev_begin, prev_end) contains all complexes of size i
  std::size_t prev_begin = 0;
  std::size_t prev_end = complexes.size();
  for (int32_t i = 1; i < l; ++i) {
    // std::cout << "Expanding complexes of size " << i << std::endl;
    for (std::size_t j = prev_begin; j < prev_end; ++j) {
      // std::cout << "Expanding complex " << complexes[j] << std::endl;
      for (const vertex_t& v : adjacent_vertices_of_complex<N>(complexes[j])) {
        // std::cout << "\nTrying adjacent vertex " << v << std::endl;
        complex_t<N> maybe_complex(complexes[j]);
        maybe_complex[v] = true;
        maybe_complex = unique_complex<N>(maybe_complex, transformations);
        // std::cout << "Maybe new complex: " << maybe_complex << std::endl;
        if (std::find(complexes.begin() + prev_end, complexes.end(),
                      maybe_complex) == complexes.end()) {
          if (is_complex<N>(maybe_complex)) {
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

template std::vector<vertex_t> adjacent_vertices_of_complex<2>(
    const complex_t<2>& complex);
template std::vector<vertex_t> adjacent_vertices_of_complex<3>(
    const complex_t<3>& complex);
template std::vector<vertex_t> adjacent_vertices_of_complex<4>(
    const complex_t<4>& complex);
template std::vector<vertex_t> adjacent_vertices_of_complex<5>(
    const complex_t<5>& complex);
template std::vector<vertex_t> adjacent_vertices_of_complex<6>(
    const complex_t<6>& complex);
template std::vector<vertex_t> adjacent_vertices_of_complex<7>(
    const complex_t<7>& complex);

template complex_t<2> unique_complex<2>(
    const complex_t<2>& complex,
    const std::vector<vertex_trans_t>& transformations);
template complex_t<3> unique_complex<3>(
    const complex_t<3>& complex,
    const std::vector<vertex_trans_t>& transformations);
template complex_t<4> unique_complex<4>(
    const complex_t<4>& complex,
    const std::vector<vertex_trans_t>& transformations);
template complex_t<5> unique_complex<5>(
    const complex_t<5>& complex,
    const std::vector<vertex_trans_t>& transformations);
template complex_t<6> unique_complex<6>(
    const complex_t<6>& complex,
    const std::vector<vertex_trans_t>& transformations);
template complex_t<7> unique_complex<7>(
    const complex_t<7>& complex,
    const std::vector<vertex_trans_t>& transformations);

template sliceable_set_t<2> complex_to_sliceable_set<2>(
    const complex_t<2>& complex, const std::vector<edge_t>& edges);
template sliceable_set_t<3> complex_to_sliceable_set<3>(
    const complex_t<3>& complex, const std::vector<edge_t>& edges);
template sliceable_set_t<4> complex_to_sliceable_set<4>(
    const complex_t<4>& complex, const std::vector<edge_t>& edges);
template sliceable_set_t<5> complex_to_sliceable_set<5>(
    const complex_t<5>& complex, const std::vector<edge_t>& edges);
template sliceable_set_t<6> complex_to_sliceable_set<6>(
    const complex_t<6>& complex, const std::vector<edge_t>& edges);
template sliceable_set_t<7> complex_to_sliceable_set<7>(
    const complex_t<7>& complex, const std::vector<edge_t>& edges);

template std::vector<sliceable_set_t<2>> complexes_to_usr<2>(
    const std::vector<complex_t<2>>& complexes,
    const std::vector<edge_trans_t>& edge_transformations,
    const std::vector<edge_t>& edges);
template std::vector<sliceable_set_t<3>> complexes_to_usr<3>(
    const std::vector<complex_t<3>>& complexes,
    const std::vector<edge_trans_t>& edge_transformations,
    const std::vector<edge_t>& edges);
template std::vector<sliceable_set_t<4>> complexes_to_usr<4>(
    const std::vector<complex_t<4>>& complexes,
    const std::vector<edge_trans_t>& edge_transformations,
    const std::vector<edge_t>& edges);
template std::vector<sliceable_set_t<5>> complexes_to_usr<5>(
    const std::vector<complex_t<5>>& complexes,
    const std::vector<edge_trans_t>& edge_transformations,
    const std::vector<edge_t>& edges);
template std::vector<sliceable_set_t<6>> complexes_to_usr<6>(
    const std::vector<complex_t<6>>& complexes,
    const std::vector<edge_trans_t>& edge_transformations,
    const std::vector<edge_t>& edges);
template std::vector<sliceable_set_t<7>> complexes_to_usr<7>(
    const std::vector<complex_t<7>>& complexes,
    const std::vector<edge_trans_t>& edge_transformations,
    const std::vector<edge_t>& edges);

template std::vector<sliceable_set_t<2>> complexes_to_mss<2>(
    const std::vector<complex_t<2>>& complexes,
    const std::vector<vertex_trans_t>& vertex_transformations,
    const std::vector<edge_t>& edges);
template std::vector<sliceable_set_t<3>> complexes_to_mss<3>(
    const std::vector<complex_t<3>>& complexes,
    const std::vector<vertex_trans_t>& vertex_transformations,
    const std::vector<edge_t>& edges);
template std::vector<sliceable_set_t<4>> complexes_to_mss<4>(
    const std::vector<complex_t<4>>& complexes,
    const std::vector<vertex_trans_t>& vertex_transformations,
    const std::vector<edge_t>& edges);
template std::vector<sliceable_set_t<5>> complexes_to_mss<5>(
    const std::vector<complex_t<5>>& complexes,
    const std::vector<vertex_trans_t>& vertex_transformations,
    const std::vector<edge_t>& edges);
template std::vector<sliceable_set_t<6>> complexes_to_mss<6>(
    const std::vector<complex_t<6>>& complexes,
    const std::vector<vertex_trans_t>& vertex_transformations,
    const std::vector<edge_t>& edges);
template std::vector<sliceable_set_t<7>> complexes_to_mss<7>(
    const std::vector<complex_t<7>>& complexes,
    const std::vector<vertex_trans_t>& vertex_transformations,
    const std::vector<edge_t>& edges);

template std::vector<complex_t<2>> compute_cut_complexes<2>(
    const std::vector<vertex_trans_t>& transformations);
template std::vector<complex_t<3>> compute_cut_complexes<3>(
    const std::vector<vertex_trans_t>& transformations);
template std::vector<complex_t<4>> compute_cut_complexes<4>(
    const std::vector<vertex_trans_t>& transformations);
template std::vector<complex_t<5>> compute_cut_complexes<5>(
    const std::vector<vertex_trans_t>& transformations);
template std::vector<complex_t<6>> compute_cut_complexes<6>(
    const std::vector<vertex_trans_t>& transformations);
template std::vector<complex_t<7>> compute_cut_complexes<7>(
    const std::vector<vertex_trans_t>& transformations);

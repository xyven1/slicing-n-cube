#ifndef N_CUBE_COMPLEX_H_
#define N_CUBE_COMPLEX_H_

#include <CGAL/Gmpz.h>
#include <CGAL/QP_functions.h>
#include <CGAL/QP_models.h>

#include <algorithm>
#include <array>
#include <bitset>
#include <cstdint>
#include <functional>
#include <vector>

#include "bitset_comparator.hpp"
#include "edge.hpp"
#include "sliceable_set.hpp"
#include "vertex.hpp"

namespace ncube {

/* complex[v] is true if vertex v is in the complex and false otherwise. */
template <int32_t N>
using complex_t = std::bitset<num_vertices(N)>;

/**
 *  Returns true if a degree one polynomial of N variables separates the given
 *  set of vertices from its complement and returns false otherwise.
 **/
template <int32_t N>
bool is_complex_degree_one(const complex_t<N>& complex) {
  using IT = int32_t;
  using ET = CGAL::Gmpz;
  using Program = CGAL::Quadratic_program<IT>;
  using Solution = CGAL::Quadratic_program_solution<ET>;
  Program lp(CGAL::SMALLER, false, 0, false, 0);
  for (vertex_t v = 0; v < num_vertices(N); ++v) {
    // The variable sign serves two purposes:
    // 1. Invert the inequality operator for vertices not part of the complex by
    //    multiplying the inequality with -1.
    // 2. Express the strict inequality `a * x < 0` as `a * x <= -eps`.
    //    To retain integer inputs the inequality is then rescaled by eps^-1.
    const int32_t sign = (complex[v]) ? 1000 : -1000;
    for (int32_t i = 0; i < N; ++i) {
      const auto coordinate_i = get_coordinate(v, i);
      lp.set_a(i, v, sign * coordinate_i);
    }
    lp.set_a(N, v, sign);
    lp.set_b(v, -1);
  }
  Solution s = CGAL::solve_linear_program(lp, ET());
  return !s.is_infeasible();
}

/**
 *  Returns true if a degree two polynomial of N variables separates the given
 *  set of vertices from its complement and returns false otherwise.
 **/
template <int32_t N>
bool is_complex_degree_two(const complex_t<N>& complex) {
  using IT = int32_t;
  using ET = CGAL::Gmpz;
  using Program = CGAL::Quadratic_program<IT>;
  using Solution = CGAL::Quadratic_program_solution<ET>;
  Program lp(CGAL::SMALLER, false, 0, false, 0);
  for (vertex_t v = 0; v < num_vertices(N); ++v) {
    // The variable sign serves two purposes:
    // 1. Invert the inequality operator for vertices not part of the complex by
    //    multiplying the inequality with -1.
    // 2. Express the strict inequality `a * x < 0` as `a * x <= -eps`.
    //    To retain integer inputs the inequality is then rescaled by eps^-1.
    const int32_t sign = (complex[v]) ? 1000 : -1000;
    for (int32_t i = 0; i < N; ++i) {
      const auto coordinate_i = get_coordinate(v, i);
      lp.set_a(i, v, sign * coordinate_i);
    }
    for (int32_t i = 0; i < N; ++i) {
      const auto coordinate_i = get_coordinate(v, i);
      for (int32_t j = 0; j < N; ++j) {
        const auto coordinate_j = get_coordinate(v, j);
        lp.set_a(N + i * N + j, v, sign * coordinate_i * coordinate_j);
      }
    }
    lp.set_a(N + N * N, v, sign);
    lp.set_b(v, -1);
  }
  Solution s = CGAL::solve_linear_program(lp, ET());
  return !s.is_infeasible();
}

/**
 *  Returns the unique symmetric representation of a cut complex.
 **/
template <int32_t N>
complex_t<N> unique_complex(const complex_t<N>& complex) {
  std::array<int32_t, N> permutation;
  for (int32_t i = 0; i < N; ++i) {
    permutation[i] = i;
  }
  // Since the unique symmetric representation of a cut complex is defined as
  // the lexicographically smallest transformation, a transformation may be
  // aborted as soon as any resulting bit (starting from the leftmost bit) is 1
  // and the corresponding bit of the current minimum is 0.
  complex_t<N> min_complex(complex);
  do {
    for (int32_t signs = 0; signs < num_vertices(N); ++signs) {
      complex_t<N> complex_trans;
      bool is_new_min = false;
      for (vertex_t v = num_vertices(N) - 1; v >= 0; --v) {
        const vertex_t v_trans = transform_vertex<N>(v, permutation, signs);
        // This actually computes the inverse transformation, but because the
        // algorithm goes through all transformations it ultimately ends up
        // being the same.
        complex_trans[v] = complex[v_trans];
        if (!is_new_min && complex_trans[v] && !min_complex[v]) {
          break;
        }
        is_new_min |= !complex_trans[v] && min_complex[v];
      }
      if (is_new_min) {
        min_complex = complex_trans;
      }
    }
  } while (std::next_permutation(permutation.begin(), permutation.end()));
  return min_complex;
}

/**
 *  Returns all vertices that are adjacent to a vertex in the cut complex but
 *  are not part of the complex already.
 **/
template <int32_t N>
std::vector<vertex_t> adjacent_vertices_of_complex(
    const complex_t<N>& complex) {
  std::vector<vertex_t> adjacent_vertices;
  for (vertex_t v = 0; v < num_vertices(N); ++v) {
    if (complex[v]) {
      for (int32_t i = 0; i < N; ++i) {
        const vertex_t neighbour = get_neighbour(v, i);
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

/**
 *  Returns the slicing corresponding to the exterior edges of a cut complex.
 **/
template <int32_t N>
sliceable_set_t<N> complex_to_sliceable_set(const complex_t<N>& complex,
                                            const edge_lexicon_t<N>& edges) {
  sliceable_set_t<N> sliceable_set;
  for (vertex_t v = 0; v < num_vertices(N); ++v) {
    if (complex[v]) {
      for (int32_t i = 0; i < N; ++i) {
        const vertex_t u = get_neighbour(v, i);
        if (!complex[u]) {
          const edge_t e = (u < v) ? edge_t(u, v) : edge_t(v, u);
          sliceable_set[edge_to_int<N>(e, edges)] = true;
        }
      }
    }
  }
  return sliceable_set;
}

/**
 *  Returns the unique symmetric representation of all slicings induced by the
 *  cut complexes.
 *
 *  The returned slicings are sorted in lexicographic order.
 **/
template <int32_t N>
std::vector<sliceable_set_t<N>> complexes_to_usr(
    const std::vector<complex_t<N>>& complexes,
    const edge_lexicon_t<N>& edges) {
  std::vector<sliceable_set_t<N>> usr;
  for (const auto& complex : complexes) {
    auto ss = complex_to_sliceable_set<N>(complex, edges);
    ss = unique_sliceable_set<N>(ss, edges);
    usr.push_back(ss);
  }
  std::sort(usr.begin(), usr.end());
  usr.erase(std::unique(usr.begin(), usr.end()), usr.end());
  return usr;
}

/**
 *  Returns the unique symmetric representation of all cut complexes subject
 *  to a given function that decides if a set of vertices is a cut complex.
 **/
template <int32_t N>
std::vector<complex_t<N>> compute_cut_complexes(
    std::function<bool(const complex_t<N>&)> is_complex) {
  // There is exactly one USR of all complexes of size 1.
  std::vector<complex_t<N>> complexes = {{1}};
  // The range [prev_begin, prev_end) contains all complexes of size i.
  std::size_t prev_begin = 0;
  std::size_t prev_end = complexes.size();
  for (int32_t i = 1; i < num_vertices(N) / 2; ++i) {
    for (std::size_t j = prev_begin; j < prev_end; ++j) {
      for (const auto& v : adjacent_vertices_of_complex<N>(complexes[j])) {
        complex_t<N> new_complex = complexes[j];
        new_complex[v] = true;
        new_complex = unique_complex<N>(new_complex);
        if (std::find(complexes.begin() + prev_end, complexes.end(),
                      new_complex) == complexes.end()) {
          if (is_complex(new_complex)) {
            complexes.push_back(new_complex);
          }
        }
      }
    }
    prev_begin = prev_end;
    prev_end = complexes.size();
  }
  return complexes;
}

/**
 *  Returns the unique symmetric representation of all cut complexes subject
 *  to a degree one polynomial (i.e. a normal hyperplane).
 **/
template <int32_t N>
std::vector<complex_t<N>> compute_cut_complexes_degree_one() {
  return compute_cut_complexes<N>(is_complex_degree_one<N>);
}

/**
 *  Returns the unique symmetric representation of all cut complexes subject
 *  to a degree two polynomial.
 **/
template <int32_t N>
std::vector<complex_t<N>> compute_cut_complexes_degree_two() {
  return compute_cut_complexes<N>(is_complex_degree_two<N>);
}

}  // namespace ncube

#endif  // N_CUBE_COMPLEX_H_

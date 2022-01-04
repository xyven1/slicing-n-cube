#ifndef N_CUBE_LP_H_
#define N_CUBE_LP_H_

#include <CGAL/Gmpz.h>
#include <CGAL/QP_functions.h>
#include <CGAL/QP_models.h>

#include <cstdint>

#include "common.hpp"

typedef int32_t IT;
typedef CGAL::Gmpz ET;

typedef CGAL::Quadratic_program<IT> Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

/**
 *  Returns true if a degree one polynomial of N variables separates the given
 *  set of vertices from its complement and returns false otherwise.
 **/
template <int32_t N>
bool is_complex_degree_one(const complex_t<N>& complex) {
  Program lp(CGAL::SMALLER, false, 0, false, 0);
  for (vertex_t v = 0; v < num_vertices(N); ++v) {
    // invert inequality operator for vertices not part of the complex
    const int32_t sign = (complex[v]) ? 1000 : -1000;
    for (int32_t i = 0; i < N; ++i) {
      int32_t val_i = (v >> i) & 1;
      int32_t coordinate_i = (val_i) ? 1 : -1;
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
  Program lp(CGAL::SMALLER, false, 0, false, 0);
  for (vertex_t v = 0; v < num_vertices(N); ++v) {
    // invert inequality operator for vertices not part of the complex
    const int32_t sign = (complex[v]) ? 1000 : -1000;
    for (int32_t i = 0; i < N; ++i) {
      const int32_t coordinate_i = ((v >> i) & 1) ? 1 : -1;
      lp.set_a(i, v, sign * coordinate_i);
    }
    for (int32_t i = 0; i < N; ++i) {
      const int32_t coordinate_i = ((v >> i) & 1) ? 1 : -1;
      for (int32_t j = 0; j < N; ++j) {
        const int32_t coordinate_j = ((v >> j) & 1) ? 1 : -1;
        lp.set_a(N + i * N + j, v, sign * coordinate_i * coordinate_j);
      }
    }
    lp.set_a(N + N * N, v, sign);
    lp.set_b(v, -1);
  }
  Solution s = CGAL::solve_linear_program(lp, ET());
  return !s.is_infeasible();
}

#endif  // N_CUBE_LP_H_

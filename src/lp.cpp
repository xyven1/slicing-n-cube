#include "lp.hpp"

#include <CGAL/Gmpz.h>
#include <CGAL/QP_functions.h>
#include <CGAL/QP_models.h>

#include <cstdint>

#include "common.hpp"

typedef int32_t IT;
typedef CGAL::Gmpz ET;

typedef CGAL::Quadratic_program<IT> Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

template <int32_t N>
bool is_complex(const complex_t<N>& complex) {
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

template bool is_complex<2>(const complex_t<2>& complex);
template bool is_complex<3>(const complex_t<3>& complex);
template bool is_complex<4>(const complex_t<4>& complex);
template bool is_complex<5>(const complex_t<5>& complex);
template bool is_complex<6>(const complex_t<6>& complex);
template bool is_complex<7>(const complex_t<7>& complex);

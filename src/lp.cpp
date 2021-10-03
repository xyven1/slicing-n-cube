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

bool is_complex(const complex_t& complex, int32_t n) {
  Program lp(CGAL::SMALLER, false, 0, false, 0);
  for (vertex_t v = 0; v < (1u << n); ++v) {
    // invert inequality operator for vertices not part of the complex
    const int32_t sign = (complex[v]) ? 1000 : -1000;
    for (int32_t i = 0; i < n; ++i) {
      uint32_t val_i = (v >> i) & 1;
      int32_t coordinate_i = (val_i) ? 1 : -1;
      lp.set_a(i, v, sign * coordinate_i);
    }
    lp.set_a(n, v, sign);
    lp.set_b(v, -1);
  }
  lp.set_l(n, true, 0);
  for (int32_t i = 0; i < n; ++i) {
    lp.set_l(i, true, 1);
    Solution s = CGAL::solve_linear_program(lp, ET());
    // std::cout << s << std::endl;
    if (!s.is_infeasible()) {
      // std::vector<CGAL::Quotient<ET>> variable_values;
      // variable_values.reserve(n + 1);
      // for (auto it = s.variable_values_begin(); it !=
      // s.variable_values_end();
      //      ++it) {
      //   variable_values.push_back(*it);
      // }
      return true;
    }
    lp.set_l(i, false);
    lp.set_u(i, true, -1);
    s = CGAL::solve_linear_program(lp, ET());
    // std::cout << s << std::endl;
    if (!s.is_infeasible()) {
      return true;
    }
    lp.set_u(i, false);
  }
  return false;
}

/*
bool verify_complex(const std::vector<CGAL::Quotient<ET>>& variable_values,
                    const complex_t& complex, int32_t n) {
  const std::vector<edge_t> edges = complex_to_edges(complex, n);
  const CGAL::Quotient<ET> zero(0);
  const CGAL::Quotient<ET> eps(1, 1000);
  for (const auto& [u, v] : edges) {
    CGAL::Quotient<ET> eval_u = variable_values[n];
    CGAL::Quotient<ET> eval_v = variable_values[n];
    for (int32_t i = 0; i < n; ++i) {
      uint32_t val_u_i = (u >> i) & 1;
      uint32_t val_v_i = (v >> i) & 1;
      int32_t coordinate_u_i = (val_u_i) ? 1 : -1;
      int32_t coordinate_v_i = (val_v_i) ? 1 : -1;
      eval_u += coordinate_u_i * variable_values[i];
      eval_v += coordinate_v_i * variable_values[i];
    }
    // we require eval_u < 0 which is emulated by eval_u + eps <= 0
    const int32_t sign_u = (complex.count(u) == 1) ? 1 : -1;
    const int32_t sign_v = (complex.count(v) == 1) ? 1 : -1;
    if (eval_u + sign_u * eps > zero) {
      return true;
    }
    if (eval_v + sign_v * eps > zero) {
      return true;
    }
  }
  return false;
}
*/

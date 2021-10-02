#include <CGAL/Gmpz.h>
#include <CGAL/QP_functions.h>
#include <CGAL/QP_models.h>

#include <algorithm>
#include <bitset>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include "prettyprint.hpp"

// typedef CGAL::Gmpq IT;
// typedef CGAL::Gmpq ET;
typedef int32_t IT;
typedef CGAL::Gmpz ET;

typedef CGAL::Quadratic_program<IT> Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

// The i-th least significant bit stores the i-th coordinate
using vertex_t = uint32_t;
using edge_t = std::pair<vertex_t, vertex_t>;
using symmetry_t = std::vector<std::pair<int32_t, int32_t>>;
using complex_t = std::set<vertex_t>;

template <typename T>
T pow(T base, T exponent) {
  T power = 1;
  for (T i = 0; i < exponent; ++i) {
    power *= base;
  }
  return power;
}

template <typename T>
T factorial(T n) {
  T f = 1;
  for (T i = 2; i <= n; ++i) {
    f *= i;
  }
  return f;
}

std::string vertex_to_str(vertex_t v, int32_t n) {
  std::string str;
  str.reserve(n);
  for (int32_t i = 0; i < n; ++i) {
    const uint32_t bit_i = (v >> i) & 1;
    str += (bit_i) ? "1" : "0";
  }
  return str;
}

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

complex_t transform_complex(const complex_t& complex, const symmetry_t& sym,
                            int32_t n) {
  complex_t transformation;
  for (const vertex_t& v : complex) {
    transformation.insert(transform_vertex(sym, v, n));
  }
  return transformation;
}

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

bool is_complex(const complex_t& complex, int32_t n, int32_t t,
                int32_t nonzero_i, bool nonzero_positive) {
  std::cout << "Verifying" << std::endl;
  constexpr int32_t scale = 1000;
  Program lp(CGAL::SMALLER, false, 0, false, 0);
  for (vertex_t v = 0; v < (1u << n); ++v) {
    // invert inequality operator for vertices not part of the complex
    const int32_t sign = (complex.count(v) == 1) ? 1 : -1;
    for (int32_t i = 0; i < n; ++i) {
      uint32_t val_i = (v >> i) & 1;
      int32_t coordinate_i = (val_i) ? 1 : -1;
      lp.set_a(i, v, sign * scale * coordinate_i);
    }
    lp.set_b(v, sign * scale * -t - 1);
  }
  if (nonzero_positive) {
    lp.set_l(nonzero_i, true, 1);
  } else {
    lp.set_u(nonzero_i, true, -1);
  }
  Solution s = CGAL::solve_linear_program(lp, ET());
  std::cout << s;
  if (!s.is_infeasible()) {
    std::cout << "  t: " << CGAL::Quotient<ET>(scale * t - 1, scale)
              << std::endl;
  }
  std::cout << std::endl;
  return !s.is_infeasible();
}

bool is_complex(const complex_t& complex, int32_t n) {
  Program lp(CGAL::SMALLER, false, 0, false, 0);
  for (vertex_t v = 0; v < (1u << n); ++v) {
    // invert inequality operator for vertices not part of the complex
    const int32_t sign = (complex.count(v) == 1) ? 1000 : -1000;
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
      // if (!detect_bullshit(variable_values, complex, n)) {
      //   return true;
      // }
      return true;
      // const CGAL::Quotient<ET> t = *(s.variable_values_begin() + n);
      // if (t.denominator() != 1) {
      //   throw std::runtime_error("Denominator is not 1");
      // }
      // const int32_t tt = static_cast<int32_t>(t.numerator().to_double());
      // if (is_complex(complex, n, tt, i, true)) {
      //   return true;
      // }
    }
    lp.set_l(i, false);
    lp.set_u(i, true, -1);
    s = CGAL::solve_linear_program(lp, ET());
    // std::cout << s << std::endl;
    if (!s.is_infeasible()) {
      // std::vector<CGAL::Quotient<ET>> variable_values;
      // variable_values.reserve(n + 1);
      // for (auto it = s.variable_values_begin(); it !=
      // s.variable_values_end();
      //      ++it) {
      //   variable_values.push_back(*it);
      // }
      // if (!detect_bullshit(variable_values, complex, n)) {
      //   return true;
      // }
      return true;
      // const CGAL::Quotient<ET> t = *(s.variable_values_begin() + n);
      // if (t.denominator() != 1) {
      //   throw std::runtime_error("Denominator is not 1");
      // }
      // const int32_t tt = static_cast<int32_t>(t.numerator().to_double());
      // if (is_complex(complex, n, tt, i, false)) {
      //   return true;
      // }
    }
    lp.set_u(i, false);
  }
  return false;
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
  const std::vector<complex_t> complexes = compute_cut_complexes(n);
  std::cout << complexes.size() << std::endl;
  // for (const complex_t& c : complexes) {
  //   std::cout << complex_to_edges(c, n) << std::endl;
  // }
}
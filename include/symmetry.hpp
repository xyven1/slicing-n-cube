#ifndef SYMMETRY_H
#define SYMMETRY_H

#include <cstdint>
#include <vector>

#include "common.hpp"

std::vector<symmetry_t> compute_symmetries(int32_t n);

vertex_t transform_vertex(const symmetry_t& sym, vertex_t v, int32_t n);

complex_t transform_complex(const complex_t& complex, const symmetry_t& sym,
                            int32_t n);

#endif

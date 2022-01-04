#ifndef N_CUBE_VERTEX_H_
#define N_CUBE_VERTEX_H_

#include <cstdint>

// The i-th least significant bit stores the i-th coordinate
using vertex_t = int32_t;

// 2^n
constexpr int32_t num_vertices(int32_t n) { return 1 << n; }

#endif  // N_CUBE_VERTEX_H_

#ifndef N_CUBE_VERTEX_H_
#define N_CUBE_VERTEX_H_

#include <array>
#include <cstdint>

namespace ncube {

/**
 *  The i-th least significant bit stores the i-th coordinate.
 *
 *  If the bit is set, the coordinate is 1. Else, the coordinate is -1.
 **/
using vertex_t = int32_t;

// 2^n
constexpr int32_t num_vertices(int32_t n) { return 1 << n; }

/**
 *  Returns the i-th coordinate of a vertex.
 **/
int32_t get_coordinate(vertex_t v, int32_t i) {
  const int32_t mask = 1 << i;
  return (v & mask) ? 1 : -1;
}

/**
 *  Returns the i-th neighbour of a vertex.
 *
 *  The i-th neighbour of a vertex is the vertex which differs in coordinate i.
 **/
vertex_t get_neighbour(vertex_t v, int32_t i) {
  const int32_t mask = 1 << i;
  return v ^ mask;
}

/**
 *  Returns a symmetric transformation of a vertex.
 *
 *  Every coordinate i of the vertex is permuted to position permutation[i] and
 *  its sign is flipped if the i-th least significant bit of signs is set.
 **/
template <int32_t N>
vertex_t transform_vertex(vertex_t v, const std::array<int32_t, N>& permutation,
                          int32_t signs) {
  vertex_t v_trans = 0;
  for (int32_t i = 0; i < N; ++i) {
    const bool change_sign = (signs >> i) & 1;
    const bool coordinate_bit = (v >> i) & 1;
    // != flips the coordinate_bit if and only if change_sign is true
    const bool new_coordinate_bit = coordinate_bit != change_sign;
    v_trans |= new_coordinate_bit << permutation[i];
  }
  return v_trans;
}

}  // namespace ncube

#endif  // N_CUBE_VERTEX_H_

#ifndef N_CUBE_LP_H_
#define N_CUBE_LP_H_

#include <cstdint>

#include "common.hpp"

/**
 *  Returns true if a degree one polynomial of N variables separates the given
 *  set of vertices from its complement and returns false otherwise.
 **/
template <int32_t N>
bool is_complex_degree_one(const complex_t<N>& complex);

/**
 *  Returns true if a degree two polynomial of N variables separates the given
 *  set of vertices from its complement and returns false otherwise.
 **/
template <int32_t N>
bool is_complex_degree_two(const complex_t<N>& complex);

#endif  // N_CUBE_LP_H_

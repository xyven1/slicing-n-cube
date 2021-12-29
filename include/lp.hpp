#ifndef LP_H
#define LP_H

#include <cstdint>

#include "common.hpp"

template <int32_t N>
bool is_complex(const complex_t<N>& complex);

template <int32_t N>
bool is_complex_degree_two(const complex_t<N>& complex);

#endif

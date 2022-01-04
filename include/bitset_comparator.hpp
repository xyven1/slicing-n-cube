#ifndef N_CUBE_BITSET_COMPARATOR_H_
#define N_CUBE_BITSET_COMPARATOR_H_

#include <bitset>

namespace std {

template <std::size_t N>
bool operator<(const std::bitset<N>& x, const std::bitset<N>& y) {
  for (std::size_t i = N - 1; i < N; --i) {
    if (x[i] != y[i]) {
      return y[i];
    }
  }
  return false;
}

template <std::size_t N>
bool operator>(const std::bitset<N>& x, const std::bitset<N>& y) {
  for (std::size_t i = N - 1; i < N; --i) {
    if (x[i] != y[i]) {
      return x[i];
    }
  }
  return false;
}

}  // namespace std

#endif  // N_CUBE_BITSET_COMPARATOR_H_

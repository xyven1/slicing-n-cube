#ifndef N_CUBE_SLICEABLE_SET_H_
#define N_CUBE_SLICEABLE_SET_H_

#include <algorithm>
#include <array>
#include <bitset>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <thread>
#include <vector>

#include "bitset_comparator.hpp"
#include "edge.hpp"
#include "symmetry.hpp"
#include "vertex.hpp"

/* edges[e] is true if edge e is in the sliceable set and false otherwise. */
template <int32_t N>
using sliceable_set_t = std::bitset<num_edges(N)>;

/**
 *  Returns a balanced distribution of m units of work across num_threads
 *  workers.
 *
 *  For example, assign_workload(15, 4) returns {4, 4, 4, 3}.
 **/
std::vector<std::size_t> assign_workload(std::size_t m,
                                         unsigned int num_threads) {
  std::vector<std::size_t> work_loads(num_threads, m / num_threads);
  for (std::size_t i = 0; i < m % num_threads; ++i) {
    ++work_loads[i];
  }
  return work_loads;
}

/**
 *  Returns the unique symmetric representation of a sliceable set.
 **/
template <int32_t N>
sliceable_set_t<N> unique_sliceable_set(const sliceable_set_t<N>& ss,
                                        const std::vector<edge_t>& edges) {
  std::array<int32_t, N> permutation;
  for (int32_t i = 0; i < N; ++i) {
    permutation[i] = i;
  }
  // Since the unique symmetric representation of a sliceable set is defined as
  // the lexicographically smallest transformation, a transformation may be
  // aborted as soon as any resulting bit (starting from the leftmost bit) is 1
  // and the corresponding bit of the current minimum is 0.
  sliceable_set_t<N> min_ss = ss;
  do {
    for (int32_t signs = 0; signs < num_vertices(N); ++signs) {
      sliceable_set_t<N> ss_trans;
      bool is_new_min = false;
      for (int32_t e = num_edges(N) - 1; e >= 0; --e) {
        const auto edge_trans = transform_edge<N>(edges[e], permutation, signs);
        const auto e_trans = edge_to_int(edge_trans, edges);
        // This actually computes the inverse transformation, but because the
        // algorithm goes through all transformations it ultimately ends up
        // being the same.
        ss_trans[e] = ss[e_trans];
        if (!is_new_min && ss_trans[e] && !min_ss[e]) {
          break;
        }
        is_new_min |= !ss_trans[e] && min_ss[e];
      }
      if (is_new_min) {
        min_ss = ss_trans;
      }
    }
  } while (std::next_permutation(permutation.begin(), permutation.end()));
  return min_ss;
}

/**
 *  Returns the unique symmetric representation of the pairwise unions of two
 *  lists of sliceable sets.
 **/
template <int32_t N>
std::vector<sliceable_set_t<N>> pairwise_unions(
    const std::vector<sliceable_set_t<N>>& sets_1,
    const std::vector<sliceable_set_t<N>>& sets_2,
    const std::vector<edge_t>& edges) {
  std::vector<sliceable_set_t<N>> unions;
  for (const auto& set_1 : sets_1) {
    for (const auto& set_2 : sets_2) {
      const auto usr = unique_sliceable_set<N>(set_1 | set_2, edges);
      const auto is_superset_of_usr = [usr](const sliceable_set_t<N>& ss) {
        return (usr | ss) == ss;
      };
      if (std::none_of(unions.begin(), unions.end(), is_superset_of_usr)) {
        const auto is_subset_of_usr = [usr](const sliceable_set_t<N>& ss) {
          return (usr | ss) == usr;
        };
        const auto it =
            remove_if(unions.begin(), unions.end(), is_subset_of_usr);
        unions.erase(it, unions.end());
        unions.push_back(usr);
      }
    }
  }
  return unions;
}

/**
 *  Stores the unique symmetric representation of the pairwise unions of two
 *  lists of sliceable sets in a range. Does NOT discard duplicates.
 **/
template <int32_t N>
void pairwise_unions_all(
    typename std::vector<sliceable_set_t<N>>::const_iterator sets_1_begin,
    typename std::vector<sliceable_set_t<N>>::const_iterator sets_1_end,
    typename std::vector<sliceable_set_t<N>>::const_iterator sets_2_begin,
    typename std::vector<sliceable_set_t<N>>::const_iterator sets_2_end,
    typename std::vector<sliceable_set_t<N>>::iterator unions_begin,
    const std::vector<edge_t>& edges) {
  auto unions_it = unions_begin;
  for (auto set_1 = sets_1_begin; set_1 != sets_1_end; ++set_1) {
    for (auto set_2 = sets_2_begin; set_2 != sets_2_end; ++set_2) {
      *unions_it = unique_sliceable_set<N>(*set_1 | *set_2, edges);
      ++unions_it;
    }
  }
}

/**
 *  Returns the unique symmetric representation of the pairwise unions of two
 *  lists of sliceable sets.
 *
 *  This function is parallelized but requires a significant amount of memory.
 **/
template <int32_t N>
std::vector<sliceable_set_t<N>> pairwise_unions_parallel(
    const std::vector<sliceable_set_t<N>>& sets_1,
    const std::vector<sliceable_set_t<N>>& sets_2,
    const std::vector<edge_t>& edges) {
  // ensure effective parallelization
  if (sets_1.size() > sets_2.size()) {
    return pairwise_unions_parallel<N>(sets_2, sets_1, edges);
  }
  // divide mss
  const unsigned int num_threads = std::thread::hardware_concurrency();
  const auto set_1_workload = assign_workload(sets_2.size(), num_threads);
  // compute the unique symmetric representation of all pairwise unions
  std::vector<sliceable_set_t<N>> unions(sets_1.size() * sets_2.size());
  std::vector<std::thread> threads;
  std::size_t prev_workload = 0;
  for (unsigned int i = 0; i < num_threads; ++i) {
    const auto unions_begin = unions.begin() + prev_workload * sets_1.size();
    const auto set_1_begin = sets_2.begin() + prev_workload;
    const auto sets_1_end = sets_2.begin() + prev_workload + set_1_workload[i];
    threads.push_back(std::thread(pairwise_unions_all<N>, sets_1.begin(),
                                  sets_1.end(), set_1_begin, sets_1_end,
                                  unions_begin, edges));
    prev_workload += set_1_workload[i];
  }
  for (auto& t : threads) {
    t.join();
  }
  // discard duplicates
  std::sort(unions.begin(), unions.end());
  auto unions_end = std::unique(unions.begin(), unions.end());
  // discard subsets
  for (auto it = unions.begin(); it != unions_end;) {
    const auto is_superset_of_it = [it](const sliceable_set_t<N>& ss) {
      return (*it | ss) == ss;
    };
    if (std::any_of(unions.begin(), it, is_superset_of_it) ||
        std::any_of(it + 1, unions_end, is_superset_of_it)) {
      --unions_end;
      *it = *unions_end;
    } else {
      ++it;
    }
  }
  unions.erase(unions_end, unions.end());
  return unions;
}

/**
 *  Returns the number of leading (leftmost) zeros in a sliceable set.
 **/
template <int32_t N>
int32_t get_leading_zeros(const sliceable_set_t<N>& ss) {
  for (std::size_t i = 0; i < ss.size(); ++i) {
    const auto rev_i = ss.size() - i - 1;
    if (ss[rev_i]) {
      return static_cast<int32_t>(i);
    }
  }
  return static_cast<int32_t>(ss.size());
}

/**
 *  Returns the number of leading (leftmost) ones in a sliceable set.
 **/
template <int32_t N>
int32_t get_leading_ones(const sliceable_set_t<N>& ss) {
  for (std::size_t i = 0; i < ss.size(); ++i) {
    const auto rev_i = ss.size() - i - 1;
    if (!ss[rev_i]) {
      return static_cast<int32_t>(i);
    }
  }
  return static_cast<int32_t>(ss.size());
}

/**
 *  Returns true if any pairwise union of two lists of sliceable sets slices all
 *  edges and false otherwise.
 *
 *  The second list is required to be sorted.
 **/
template <int32_t N>
bool pairwise_unions_slice_cube(const std::vector<sliceable_set_t<N>>& usr,
                                const std::vector<sliceable_set_t<N>>& mss) {
  bool slices_all = false;
  for (const sliceable_set_t<N>& set_1 : usr) {
    const int32_t leading_zeros = get_leading_zeros<N>(set_1);
    for (auto set_2_it = mss.rbegin(); set_2_it != mss.rend(); ++set_2_it) {
      const int32_t leading_ones = get_leading_ones<N>(*set_2_it);
      if (leading_ones < leading_zeros) {
        break;
      }
      const auto set_union = set_1 | *set_2_it;
      slices_all |= set_union.all();
    }
  }
  return slices_all;
}

/**
 *  Returns the symmetry expansions of unique symmetric representations of
 *  sliceable sets.
 *
 *  The returned sliceable sets are sorted in lexicographic order.
 **/
template <int32_t N>
std::vector<sliceable_set_t<N>> usr_to_mss(
    const std::vector<sliceable_set_t<N>>& usr,
    const std::vector<edge_t>& edges) {
  std::vector<sliceable_set_t<N>> mss;
  mss.reserve(usr.size() * num_symmetries(N));
  for (const auto& ss : usr) {
    std::array<int32_t, N> permutation;
    for (int32_t i = 0; i < N; ++i) {
      permutation[i] = i;
    }
    do {
      for (int32_t signs = 0; signs < num_vertices(N); ++signs) {
        sliceable_set_t<N> ss_trans;
        for (int32_t e = 0; e < num_edges(N); ++e) {
          const auto edge_trans =
              transform_edge<N>(edges[e], permutation, signs);
          const auto e_trans = edge_to_int(edge_trans, edges);
          ss_trans[e_trans] = ss[e];
        }
        mss.push_back(ss_trans);
      }
    } while (std::next_permutation(permutation.begin(), permutation.end()));
  }
  std::sort(mss.begin(), mss.end());
  mss.erase(std::unique(mss.begin(), mss.end()), mss.end());
  return mss;
}

/**
 *  Stores the bitstring encoding of a sliceable set in a byte array.
 **/
template <int32_t N>
void sliceable_set_to_bytes(const sliceable_set_t<N>& ss, char* bytes,
                            std::size_t num_bytes) {
  for (std::size_t i = 0; i < num_bytes; ++i) {
    const std::size_t rev_i = num_bytes - i - 1;
    bytes[rev_i] = 0;
    for (std::size_t j = 0; j < 8; ++j) {
      if (ss[i * 8 + j]) {
        bytes[rev_i] |= static_cast<char>(1 << j);
      }
    }
  }
}

/**
 *  Returns a sliceable set from to its bitstring encoding in a byte array.
 **/
template <int32_t N>
sliceable_set_t<N> bytes_to_sliceable_set(char* bytes, std::size_t num_bytes) {
  sliceable_set_t<N> ss;
  for (std::size_t i = 0; i < num_bytes; ++i) {
    const std::size_t rev_i = num_bytes - i - 1;
    for (std::size_t j = 0; j < 8; ++j) {
      const std::size_t ss_index = i * 8 + j;
      if (ss_index < ss.size()) {
        ss[ss_index] = bytes[rev_i] & (1 << j);
      }
    }
  }
  return ss;
}

/**
 *  Returns the minimum number of bytes needed to represent n bits.
 **/
constexpr std::size_t min_bytes_to_represent_bits(std::size_t n) {
  if (n % 8 == 0) {
    return n / 8;
  }
  return n / 8 + 1;
}

/**
 *  Writes sliceable sets to a file at the given path.
 **/
template <int32_t N>
void write_to_file(const std::vector<sliceable_set_t<N>>& sets,
                   const std::filesystem::path& path) {
  constexpr auto num_bytes_sliceable_set =
      min_bytes_to_represent_bits(sliceable_set_t<N>().size());
  std::ofstream file(path, std::ios::binary);
  for (const sliceable_set_t<N>& set : sets) {
    char bytes[num_bytes_sliceable_set];
    sliceable_set_to_bytes<N>(set, bytes, num_bytes_sliceable_set);
    file.write(bytes, num_bytes_sliceable_set);
  }
}

/**
 *  Returns the sliceable sets stored in a file at the given path.
 *
 *  Naturally, this function should be called on a file created by
 *  `write_to_file`.
 **/
template <int32_t N>
std::vector<sliceable_set_t<N>> read_from_file(
    const std::filesystem::path& path) {
  constexpr auto num_bytes_sliceable_set =
      min_bytes_to_represent_bits(sliceable_set_t<N>().size());
  const auto filesize = std::filesystem::file_size(path);
  std::vector<sliceable_set_t<N>> sets;
  sets.reserve(filesize / num_bytes_sliceable_set);
  std::ifstream file(path, std::ios::binary);
  for (std::size_t i = 0; i < filesize / num_bytes_sliceable_set; ++i) {
    char bytes[num_bytes_sliceable_set];
    file.read(bytes, num_bytes_sliceable_set);
    sets.push_back(bytes_to_sliceable_set<N>(bytes, num_bytes_sliceable_set));
  }
  return sets;
}

#endif  // N_CUBE_SLICEABLE_SET_H_

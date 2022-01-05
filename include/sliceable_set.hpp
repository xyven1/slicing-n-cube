#ifndef N_CUBE_SLICEABLE_SET_H_
#define N_CUBE_SLICEABLE_SET_H_

#include <algorithm>
#include <array>
#include <bitset>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <future>
#include <vector>

#include "bitset_comparator.hpp"
#include "edge.hpp"
#include "symmetry.hpp"
#include "vertex.hpp"

// edges[e] is true iff edge e is in the sliceable set
template <int32_t N>
using sliceable_set_t = std::bitset<num_edges(N)>;

/**
 *  Returns a balanced distribution of n units of work across num_threads
 *  workers.
 *
 *  For example, assign_workload(15, 4) returns {4, 4, 4, 3}.
 **/
std::vector<std::size_t> assign_workload(std::size_t n,
                                         unsigned int num_threads) {
  std::vector<std::size_t> work_loads(num_threads, n / num_threads);
  for (std::size_t i = 0; i < n % num_threads; ++i) {
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
  sliceable_set_t<N> min_ss(ss);
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
          // min_ss is still smaller
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
 *  Returns the unique symmetric representation of the pairwise unions of
 *  maximal sliceable sets and their unique symmetric representation.
 **/
template <int32_t N>
std::vector<sliceable_set_t<N>> combine_usr_mss(
    const std::vector<sliceable_set_t<N>>& usr,
    const std::vector<sliceable_set_t<N>>& mss,
    const std::vector<edge_t>& edges) {
  std::vector<sliceable_set_t<N>> combos;
  for (const sliceable_set_t<N>& set_1 : usr) {
    for (const sliceable_set_t<N>& set_2 : mss) {
      sliceable_set_t<N> combo = set_1 | set_2;
      combo = unique_sliceable_set<N>(combo, edges);
      const auto is_superset = [combo](const sliceable_set_t<N>& other) {
        return (other | combo) == other;
      };
      if (std::none_of(combos.begin(), combos.end(), is_superset)) {
        const auto is_subset = [combo](const sliceable_set_t<N>& other) {
          return (other | combo) == combo;
        };
        const auto it = remove_if(combos.begin(), combos.end(), is_subset);
        combos.erase(it, combos.end());
        combos.push_back(combo);
      }
    }
  }
  return combos;
}

/**
 *  Stores the unique symmetric representation of the pairwise unions of
 *  maximal sliceable sets and their unique symmetric representation in a range.
 **/
template <int32_t N>
void combine_usr_mss_all(
    typename std::vector<sliceable_set_t<N>>::const_iterator usr_begin,
    typename std::vector<sliceable_set_t<N>>::const_iterator usr_end,
    typename std::vector<sliceable_set_t<N>>::const_iterator mss_begin,
    typename std::vector<sliceable_set_t<N>>::const_iterator mss_end,
    typename std::vector<sliceable_set_t<N>>::iterator combos_begin,
    const std::vector<edge_t>& edges) {
  auto combos_it = combos_begin;
  for (auto usr_it = usr_begin; usr_it != usr_end; ++usr_it) {
    for (auto mss_it = mss_begin; mss_it != mss_end; ++mss_it) {
      const auto combo = *usr_it | *mss_it;
      *combos_it = unique_sliceable_set<N>(combo, edges);
      ++combos_it;
    }
  }
}

/**
 *  Returns the unique symmetric representation of the pairwise unions of
 *  maximal sliceable sets and their unique symmetric representation.
 **/
template <int32_t N>
std::vector<sliceable_set_t<N>> combine_usr_mss_parallel(
    const std::vector<sliceable_set_t<N>>& usr,
    const std::vector<sliceable_set_t<N>>& mss,
    const std::vector<edge_t>& edges) {
  // ensure effective parallelization
  if (usr.size() > mss.size()) {
    return combine_usr_mss_parallel<N>(mss, usr, edges);
  }
  // divide mss
  const unsigned int num_threads = std::thread::hardware_concurrency();
  const auto mss_workload = assign_workload(mss.size(), num_threads);
  // compute the unique symmetric representation of all pairwise unions
  std::vector<sliceable_set_t<N>> combos(usr.size() * mss.size());
  std::vector<std::thread> threads;
  std::size_t mss_workload_prior = 0;
  for (unsigned int i = 0; i < num_threads; ++i) {
    const auto combos_begin = combos.begin() + mss_workload_prior * usr.size();
    const auto mss_begin = mss.begin() + mss_workload_prior;
    const auto mss_end = mss.begin() + mss_workload_prior + mss_workload[i];
    threads.push_back(std::thread(combine_usr_mss_all<N>, usr.begin(),
                                  usr.end(), mss_begin, mss_end, combos_begin,
                                  edges));
    mss_workload_prior += mss_workload[i];
  }
  for (auto& t : threads) {
    t.join();
  }
  // discard duplicates
  std::sort(combos.begin(), combos.end());
  auto combos_end = std::unique(combos.begin(), combos.end());
  // discard subsets
  for (auto it = combos.begin(); it != combos_end;) {
    const auto is_superset = [it](const sliceable_set_t<N>& ss) {
      return (*it | ss) == ss;
    };
    if (std::any_of(combos.begin(), it, is_superset) ||
        std::any_of(it + 1, combos_end, is_superset)) {
      --combos_end;
      *it = *combos_end;
    } else {
      ++it;
    }
  }
  combos.erase(combos_end, combos.end());
  return combos;
}

/**
 *  Returns the number of leading zeros in a sliceable set.
 **/
template <int32_t N>
int32_t get_leading_zeros(const sliceable_set_t<N>& ss) {
  for (std::size_t i = 0; i < ss.size(); ++i) {
    const std::size_t rev_i = ss.size() - i - 1;
    if (ss[rev_i]) {
      return static_cast<int32_t>(i);
    }
  }
  return static_cast<int32_t>(ss.size());
}

/**
 *  Returns the number of leading ones in a sliceable set.
 **/
template <int32_t N>
int32_t get_leading_ones(const sliceable_set_t<N>& ss) {
  for (std::size_t i = 0; i < ss.size(); ++i) {
    const std::size_t rev_i = ss.size() - i - 1;
    if (!ss[rev_i]) {
      return static_cast<int32_t>(i);
    }
  }
  return static_cast<int32_t>(ss.size());
}

/**
 *  Returns whether any pairwise union of maximal sliceable sets and their
 *  unique symmetric representation slices all edges.
 *
 *  The maximal sliceable sets are required to be sorted.
 **/
template <int32_t N>
bool combine_usr_mss_final(const std::vector<sliceable_set_t<N>>& usr,
                           const std::vector<sliceable_set_t<N>>& mss) {
  bool slices_all = false;
  for (const sliceable_set_t<N>& set_1 : usr) {
    const int32_t leading_zeros = get_leading_zeros<N>(set_1);
    for (auto set_2_it = mss.rbegin(); set_2_it != mss.rend(); ++set_2_it) {
      const int32_t leading_ones = get_leading_ones<N>(*set_2_it);
      if (leading_ones < leading_zeros) {
        break;
      }
      const sliceable_set_t<N> combo = set_1 | *set_2_it;
      slices_all |= combo.all();
    }
  }
  return slices_all;
}

/**
 *  Returns whether any pairwise union of maximal sliceable sets and their
 *  unique symmetric representation slices all edges.
 *
 *  This is an inefficient implementation that exist only for reference.
 **/
template <int32_t N>
bool combine_usr_mss_final_naive(const std::vector<sliceable_set_t<N>>& usr,
                                 const std::vector<sliceable_set_t<N>>& mss) {
  for (const sliceable_set_t<N>& set_1 : mss) {
    for (const sliceable_set_t<N>& set_2 : usr) {
      if ((set_1 | set_2).all()) {
        return true;
      }
    }
  }
  return false;
}

/**
 *  Returns the symmetry expansions of the unique symmetric representation of
 *  sliceable sets.
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
  constexpr std::size_t num_bytes_sliceable_set =
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
  constexpr std::size_t num_bytes_sliceable_set =
      min_bytes_to_represent_bits(sliceable_set_t<N>().size());
  const std::size_t filesize = std::filesystem::file_size(path);
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

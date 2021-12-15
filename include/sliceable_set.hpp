#ifndef SLICEABLE_SET_H
#define SLICEABLE_SET_H

#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <future>
#include <vector>

#include "common.hpp"

std::vector<std::size_t> assign_workload(std::size_t n,
                                         unsigned int num_threads) {
  std::vector<std::size_t> work_loads(num_threads, n / num_threads);
  for (std::size_t i = 0; i < n % num_threads; ++i) {
    ++work_loads[i];
  }
  return work_loads;
}

template <int32_t N>
bool is_subset(const sliceable_set_t<N>& subset,
               const std::vector<sliceable_set_t<N>>& combos) {
  for (const sliceable_set_t<N>& combo : combos) {
    if ((subset | combo) == combo) {
      return true;
    }
  }
  return false;
}

template <int32_t N>
sliceable_set_t<N> unique_sliceable_set_naive(
    const sliceable_set_t<N>& ss,
    const std::vector<edge_trans_t>& transformations) {
  sliceable_set_t<N> min_ss(ss);
  for (const edge_trans_t& transformation : transformations) {
    sliceable_set_t<N> ss_trans;
    for (int32_t e = 0; e < num_edges(N); ++e) {
      const int32_t e_inv = transformation[e];
      ss_trans[e] = ss[e_inv];
    }
    if (min_ss > ss_trans) {
      min_ss = ss_trans;
    }
  }
  return min_ss;
}

template <int32_t N>
sliceable_set_t<N> unique_sliceable_set(
    const sliceable_set_t<N>& ss,
    const std::vector<edge_trans_t>& transformations) {
  sliceable_set_t<N> min_ss(ss);
  for (const edge_trans_t& transformation : transformations) {
    sliceable_set_t<N> ss_trans;
    bool is_new_min = false;
    for (int32_t e = num_edges(N) - 1; e >= 0; --e) {
      const int32_t e_inversion = transformation[e];
      ss_trans[e] = ss[e_inversion];
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
  return min_ss;
}

template <int32_t N>
vertex_t transform_vertex_inv(vertex_t v,
                              const std::array<int32_t, N>& positions,
                              int32_t signs) {
  vertex_t transformation = 0;
  for (int32_t i = 0; i < N; ++i) {
    const int32_t sign_i = (signs >> i) & 1;
    const int32_t val_i = (v >> positions[i]) & 1;
    const int32_t val_i_sign = val_i ^ sign_i;
    const int32_t val_i_sign_and_position = val_i_sign << i;
    transformation |= val_i_sign_and_position;
  }
  return transformation;
}

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
        const vertex_t u =
            transform_vertex_inv<N>(edges[e].first, permutation, signs);
        const vertex_t v =
            transform_vertex_inv<N>(edges[e].second, permutation, signs);
        const edge_t edge_inv = (u < v) ? edge_t(u, v) : edge_t(v, u);
        ss_trans[e] = ss[edge_to_int(edge_inv, edges)];
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

template <int32_t N>
sliceable_set_t<N> unique_sliceable_set_seq(
    const sliceable_set_t<N>& ss,
    const std::vector<edge_trans_t>::const_iterator transformations_begin,
    const std::vector<edge_trans_t>::const_iterator transformations_end) {
  sliceable_set_t<N> min_ss(ss);
  for (auto it = transformations_begin; it != transformations_end; ++it) {
    sliceable_set_t<N> ss_trans;
    bool is_new_min = false;
    for (int32_t e = num_edges(N) - 1; e >= 0; --e) {
      const int32_t e_inversion = (*it)[e];
      ss_trans[e] = ss[e_inversion];
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
  return min_ss;
}

template <int32_t N>
sliceable_set_t<N> unique_sliceable_set_parallel(
    const sliceable_set_t<N>& ss,
    const std::vector<edge_trans_t>& transformations) {
  const auto num_threads = std::thread::hardware_concurrency();
  const auto work_size = transformations.size() / num_threads;
  std::vector<std::vector<edge_trans_t>::const_iterator> separators;
  for (unsigned int i = 0; i < num_threads; ++i) {
    separators.push_back(transformations.begin() + i * work_size);
  }
  separators.push_back(transformations.end());
  std::vector<std::future<sliceable_set_t<N>>> futures;
  for (unsigned int i = 0; i < num_threads; ++i) {
    futures.push_back(std::async(unique_sliceable_set_seq<N>, ss, separators[i],
                                 separators[i + 1]));
  }
  sliceable_set_t<N> min_ss(ss);
  for (auto& future : futures) {
    const auto future_ss = future.get();
    if (min_ss > future_ss) {
      min_ss = future_ss;
    }
  }
  return min_ss;
}

template <int32_t N>
std::vector<sliceable_set_t<N>> combine_usr_mss(
    const std::vector<sliceable_set_t<N>>& usr,
    const std::vector<sliceable_set_t<N>>& mss,
    const std::vector<edge_trans_t>& transformations) {
  std::vector<sliceable_set_t<N>> combos;
  for (const sliceable_set_t<N>& set_1 : usr) {
    for (const sliceable_set_t<N>& set_2 : mss) {
      sliceable_set_t<N> combo = set_1 | set_2;
      combo = unique_sliceable_set<N>(combo, transformations);
      if (!is_subset<N>(combo, combos)) {
        const auto p = [combo](const sliceable_set_t<N>& ss) {
          return (ss | combo) == combo;
        };
        const auto it = remove_if(combos.begin(), combos.end(), p);
        combos.erase(it, combos.end());
        combos.push_back(combo);
      }
    }
  }
  return combos;
}

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
      if (!is_subset<N>(combo, combos)) {
        const auto p = [combo](const sliceable_set_t<N>& ss) {
          return (ss | combo) == combo;
        };
        const auto it = remove_if(combos.begin(), combos.end(), p);
        combos.erase(it, combos.end());
        combos.push_back(combo);
      }
    }
  }
  return combos;
}

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

template <int32_t N>
std::vector<sliceable_set_t<N>> usr_to_mss(
    const std::vector<sliceable_set_t<N>>& usr,
    const std::vector<edge_trans_t>& transformations) {
  std::vector<sliceable_set_t<N>> mss;
  mss.reserve(usr.size() * num_symmetries(N));
  for (const sliceable_set_t<N>& ss : usr) {
    for (const edge_trans_t& transformation : transformations) {
      sliceable_set_t<N> ss_trans;
      for (int32_t e = 0; e < num_edges(N); ++e) {
        const int32_t e_inversion = transformation[e];
        ss_trans[e] = ss[e_inversion];
      }
      mss.push_back(ss_trans);
    }
  }
  std::sort(mss.begin(), mss.end());
  mss.erase(std::unique(mss.begin(), mss.end()), mss.end());
  // mss.shrink_to_fit();
  return mss;
}

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

#endif

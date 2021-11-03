#include "sliceable_set.hpp"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <vector>

#include "common.hpp"

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

template bool is_subset<2>(const sliceable_set_t<2>& subset,
                           const std::vector<sliceable_set_t<2>>& combos);
template bool is_subset<3>(const sliceable_set_t<3>& subset,
                           const std::vector<sliceable_set_t<3>>& combos);
template bool is_subset<4>(const sliceable_set_t<4>& subset,
                           const std::vector<sliceable_set_t<4>>& combos);
template bool is_subset<5>(const sliceable_set_t<5>& subset,
                           const std::vector<sliceable_set_t<5>>& combos);
template bool is_subset<6>(const sliceable_set_t<6>& subset,
                           const std::vector<sliceable_set_t<6>>& combos);
template bool is_subset<7>(const sliceable_set_t<7>& subset,
                           const std::vector<sliceable_set_t<7>>& combos);

template sliceable_set_t<2> unique_sliceable_set_naive<2>(
    const sliceable_set_t<2>& ss,
    const std::vector<edge_trans_t>& transformations);
template sliceable_set_t<3> unique_sliceable_set_naive<3>(
    const sliceable_set_t<3>& ss,
    const std::vector<edge_trans_t>& transformations);
template sliceable_set_t<4> unique_sliceable_set_naive<4>(
    const sliceable_set_t<4>& ss,
    const std::vector<edge_trans_t>& transformations);
template sliceable_set_t<5> unique_sliceable_set_naive<5>(
    const sliceable_set_t<5>& ss,
    const std::vector<edge_trans_t>& transformations);
template sliceable_set_t<6> unique_sliceable_set_naive<6>(
    const sliceable_set_t<6>& ss,
    const std::vector<edge_trans_t>& transformations);
template sliceable_set_t<7> unique_sliceable_set_naive<7>(
    const sliceable_set_t<7>& ss,
    const std::vector<edge_trans_t>& transformations);

template sliceable_set_t<2> unique_sliceable_set<2>(
    const sliceable_set_t<2>& ss,
    const std::vector<edge_trans_t>& transformations);
template sliceable_set_t<3> unique_sliceable_set<3>(
    const sliceable_set_t<3>& ss,
    const std::vector<edge_trans_t>& transformations);
template sliceable_set_t<4> unique_sliceable_set<4>(
    const sliceable_set_t<4>& ss,
    const std::vector<edge_trans_t>& transformations);
template sliceable_set_t<5> unique_sliceable_set<5>(
    const sliceable_set_t<5>& ss,
    const std::vector<edge_trans_t>& transformations);
template sliceable_set_t<6> unique_sliceable_set<6>(
    const sliceable_set_t<6>& ss,
    const std::vector<edge_trans_t>& transformations);
template sliceable_set_t<7> unique_sliceable_set<7>(
    const sliceable_set_t<7>& ss,
    const std::vector<edge_trans_t>& transformations);

template std::vector<sliceable_set_t<2>> combine_usr_mss<2>(
    const std::vector<sliceable_set_t<2>>& usr,
    const std::vector<sliceable_set_t<2>>& mss,
    const std::vector<edge_trans_t>& transformations);
template std::vector<sliceable_set_t<3>> combine_usr_mss<3>(
    const std::vector<sliceable_set_t<3>>& usr,
    const std::vector<sliceable_set_t<3>>& mss,
    const std::vector<edge_trans_t>& transformations);
template std::vector<sliceable_set_t<4>> combine_usr_mss<4>(
    const std::vector<sliceable_set_t<4>>& usr,
    const std::vector<sliceable_set_t<4>>& mss,
    const std::vector<edge_trans_t>& transformations);
template std::vector<sliceable_set_t<5>> combine_usr_mss<5>(
    const std::vector<sliceable_set_t<5>>& usr,
    const std::vector<sliceable_set_t<5>>& mss,
    const std::vector<edge_trans_t>& transformations);
template std::vector<sliceable_set_t<6>> combine_usr_mss<6>(
    const std::vector<sliceable_set_t<6>>& usr,
    const std::vector<sliceable_set_t<6>>& mss,
    const std::vector<edge_trans_t>& transformations);
template std::vector<sliceable_set_t<7>> combine_usr_mss<7>(
    const std::vector<sliceable_set_t<7>>& usr,
    const std::vector<sliceable_set_t<7>>& mss,
    const std::vector<edge_trans_t>& transformations);

template bool combine_usr_mss_final<2>(
    const std::vector<sliceable_set_t<2>>& usr,
    const std::vector<sliceable_set_t<2>>& mss);
template bool combine_usr_mss_final<3>(
    const std::vector<sliceable_set_t<3>>& usr,
    const std::vector<sliceable_set_t<3>>& mss);
template bool combine_usr_mss_final<4>(
    const std::vector<sliceable_set_t<4>>& usr,
    const std::vector<sliceable_set_t<4>>& mss);
template bool combine_usr_mss_final<5>(
    const std::vector<sliceable_set_t<5>>& usr,
    const std::vector<sliceable_set_t<5>>& mss);
template bool combine_usr_mss_final<6>(
    const std::vector<sliceable_set_t<6>>& usr,
    const std::vector<sliceable_set_t<6>>& mss);
template bool combine_usr_mss_final<7>(
    const std::vector<sliceable_set_t<7>>& usr,
    const std::vector<sliceable_set_t<7>>& mss);

template bool combine_usr_mss_final_naive<2>(
    const std::vector<sliceable_set_t<2>>& usr,
    const std::vector<sliceable_set_t<2>>& mss);
template bool combine_usr_mss_final_naive<3>(
    const std::vector<sliceable_set_t<3>>& usr,
    const std::vector<sliceable_set_t<3>>& mss);
template bool combine_usr_mss_final_naive<4>(
    const std::vector<sliceable_set_t<4>>& usr,
    const std::vector<sliceable_set_t<4>>& mss);
template bool combine_usr_mss_final_naive<5>(
    const std::vector<sliceable_set_t<5>>& usr,
    const std::vector<sliceable_set_t<5>>& mss);
template bool combine_usr_mss_final_naive<6>(
    const std::vector<sliceable_set_t<6>>& usr,
    const std::vector<sliceable_set_t<6>>& mss);
template bool combine_usr_mss_final_naive<7>(
    const std::vector<sliceable_set_t<7>>& usr,
    const std::vector<sliceable_set_t<7>>& mss);

template std::vector<sliceable_set_t<2>> usr_to_mss<2>(
    const std::vector<sliceable_set_t<2>>& usr,
    const std::vector<edge_trans_t>& transformations);
template std::vector<sliceable_set_t<3>> usr_to_mss<3>(
    const std::vector<sliceable_set_t<3>>& usr,
    const std::vector<edge_trans_t>& transformations);
template std::vector<sliceable_set_t<4>> usr_to_mss<4>(
    const std::vector<sliceable_set_t<4>>& usr,
    const std::vector<edge_trans_t>& transformations);
template std::vector<sliceable_set_t<5>> usr_to_mss<5>(
    const std::vector<sliceable_set_t<5>>& usr,
    const std::vector<edge_trans_t>& transformations);
template std::vector<sliceable_set_t<6>> usr_to_mss<6>(
    const std::vector<sliceable_set_t<6>>& usr,
    const std::vector<edge_trans_t>& transformations);
template std::vector<sliceable_set_t<7>> usr_to_mss<7>(
    const std::vector<sliceable_set_t<7>>& usr,
    const std::vector<edge_trans_t>& transformations);

template int32_t get_leading_zeros<2>(const sliceable_set_t<2>& ss);
template int32_t get_leading_zeros<3>(const sliceable_set_t<3>& ss);
template int32_t get_leading_zeros<4>(const sliceable_set_t<4>& ss);
template int32_t get_leading_zeros<5>(const sliceable_set_t<5>& ss);
template int32_t get_leading_zeros<6>(const sliceable_set_t<6>& ss);
template int32_t get_leading_zeros<7>(const sliceable_set_t<7>& ss);

template int32_t get_leading_ones<2>(const sliceable_set_t<2>& ss);
template int32_t get_leading_ones<3>(const sliceable_set_t<3>& ss);
template int32_t get_leading_ones<4>(const sliceable_set_t<4>& ss);
template int32_t get_leading_ones<5>(const sliceable_set_t<5>& ss);
template int32_t get_leading_ones<6>(const sliceable_set_t<6>& ss);
template int32_t get_leading_ones<7>(const sliceable_set_t<7>& ss);

template void write_to_file<2>(const std::vector<sliceable_set_t<2>>& sets,
                               const std::filesystem::path& path);
template void write_to_file<3>(const std::vector<sliceable_set_t<3>>& sets,
                               const std::filesystem::path& path);
template void write_to_file<4>(const std::vector<sliceable_set_t<4>>& sets,
                               const std::filesystem::path& path);
template void write_to_file<5>(const std::vector<sliceable_set_t<5>>& sets,
                               const std::filesystem::path& path);
template void write_to_file<6>(const std::vector<sliceable_set_t<6>>& sets,
                               const std::filesystem::path& path);
template void write_to_file<7>(const std::vector<sliceable_set_t<7>>& sets,
                               const std::filesystem::path& path);

template std::vector<sliceable_set_t<2>> read_from_file<2>(
    const std::filesystem::path& path);
template std::vector<sliceable_set_t<3>> read_from_file<3>(
    const std::filesystem::path& path);
template std::vector<sliceable_set_t<4>> read_from_file<4>(
    const std::filesystem::path& path);
template std::vector<sliceable_set_t<5>> read_from_file<5>(
    const std::filesystem::path& path);
template std::vector<sliceable_set_t<6>> read_from_file<6>(
    const std::filesystem::path& path);
template std::vector<sliceable_set_t<7>> read_from_file<7>(
    const std::filesystem::path& path);

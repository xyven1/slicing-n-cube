#include "sliceable_set.hpp"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <vector>

#include "common.hpp"

bool is_subset(const sliceable_set_t& subset,
               const std::vector<sliceable_set_t>& combos) {
  for (const sliceable_set_t& combo : combos) {
    if ((subset | combo) == combo) {
      return true;
    }
  }
  return false;
}

sliceable_set_t unique_sliceable_set_naive(
    const sliceable_set_t& ss, const std::vector<edge_trans_t>& transformations,
    int32_t n) {
  sliceable_set_t min_ss(ss);
  for (const edge_trans_t& transformation : transformations) {
    sliceable_set_t ss_trans;
    for (int32_t e = 0; e < num_edges(n); ++e) {
      const int32_t e_inv = transformation[e];
      ss_trans[e] = ss[e_inv];
    }
    if (min_ss > ss_trans) {
      min_ss = ss_trans;
    }
  }
  return min_ss;
}

sliceable_set_t unique_sliceable_set(
    const sliceable_set_t& ss, const std::vector<edge_trans_t>& transformations,
    int32_t n) {
  sliceable_set_t min_ss(ss);
  for (const edge_trans_t& transformation : transformations) {
    sliceable_set_t ss_trans;
    bool is_new_min = false;
    for (int32_t e = num_edges(n) - 1; e >= 0; --e) {
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

std::vector<sliceable_set_t> combine_usr_mss(
    const std::vector<sliceable_set_t>& usr,
    const std::vector<sliceable_set_t>& mss,
    const std::vector<edge_trans_t>& transformations, int32_t n) {
  std::vector<sliceable_set_t> combos;
  for (const sliceable_set_t& set_1 : usr) {
    for (const sliceable_set_t& set_2 : mss) {
      sliceable_set_t combo = set_1 | set_2;
      combo = unique_sliceable_set(combo, transformations, n);
      if (!is_subset(combo, combos)) {
        const auto p = [combo](const sliceable_set_t& ss) {
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

bool combine_usr_mss_final(const std::vector<sliceable_set_t>& usr,
                           const std::vector<sliceable_set_t>& mss) {
  bool slices_all = false;
  for (const sliceable_set_t& set_1 : usr) {
    const int32_t leading_zeros = get_leading_zeros(set_1);
    for (auto set_2_it = mss.rbegin(); set_2_it != mss.rend(); ++set_2_it) {
      const int32_t leading_ones = get_leading_ones(*set_2_it);
      if (leading_ones < leading_zeros) {
        break;
      }
      const sliceable_set_t combo = set_1 | *set_2_it;
      slices_all |= combo.all();
    }
  }
  return slices_all;
}

bool combine_usr_mss_final_naive(const std::vector<sliceable_set_t>& usr,
                                 const std::vector<sliceable_set_t>& mss) {
  for (const sliceable_set_t& set_1 : mss) {
    for (const sliceable_set_t& set_2 : usr) {
      if ((set_1 | set_2).all()) {
        return true;
      }
    }
  }
  return false;
}

std::vector<sliceable_set_t> usr_to_mss(
    const std::vector<sliceable_set_t>& usr,
    const std::vector<edge_trans_t>& transformations, int32_t n) {
  std::vector<sliceable_set_t> mss;
  mss.reserve(usr.size() * num_symmetries(n));
  for (const sliceable_set_t& ss : usr) {
    for (const edge_trans_t& transformation : transformations) {
      sliceable_set_t ss_trans;
      for (int32_t e = 0; e < num_edges(n); ++e) {
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

int32_t get_leading_zeros(const sliceable_set_t& ss) {
  for (std::size_t i = 0; i < ss.size(); ++i) {
    const std::size_t rev_i = ss.size() - i - 1;
    if (ss[rev_i]) {
      return static_cast<int32_t>(i);
    }
  }
  return static_cast<int32_t>(ss.size());
}

int32_t get_leading_ones(const sliceable_set_t& ss) {
  for (std::size_t i = 0; i < ss.size(); ++i) {
    const std::size_t rev_i = ss.size() - i - 1;
    if (!ss[rev_i]) {
      return static_cast<int32_t>(i);
    }
  }
  return static_cast<int32_t>(ss.size());
}

void sliceable_set_to_bytes(const sliceable_set_t& ss, char* bytes,
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

sliceable_set_t bytes_to_sliceable_set(char* bytes, std::size_t num_bytes) {
  sliceable_set_t ss;
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

void write_to_file(const std::vector<sliceable_set_t>& sets,
                   const std::filesystem::path& path) {
  constexpr std::size_t num_bytes_sliceable_set =
      min_bytes_to_represent_bits(sliceable_set_t().size());
  std::ofstream file(path, std::ios::binary);
  for (const sliceable_set_t& set : sets) {
    char bytes[num_bytes_sliceable_set];
    sliceable_set_to_bytes(set, bytes, num_bytes_sliceable_set);
    file.write(bytes, num_bytes_sliceable_set);
  }
}

std::vector<sliceable_set_t> read_from_file(const std::filesystem::path& path) {
  constexpr std::size_t num_bytes_sliceable_set =
      min_bytes_to_represent_bits(sliceable_set_t().size());
  const std::size_t filesize = std::filesystem::file_size(path);
  std::vector<sliceable_set_t> sets;
  sets.reserve(filesize / num_bytes_sliceable_set);
  std::ifstream file(path, std::ios::binary);
  for (std::size_t i = 0; i < filesize / num_bytes_sliceable_set; ++i) {
    char bytes[num_bytes_sliceable_set];
    file.read(bytes, num_bytes_sliceable_set);
    sets.push_back(bytes_to_sliceable_set(bytes, num_bytes_sliceable_set));
  }
  return sets;
}

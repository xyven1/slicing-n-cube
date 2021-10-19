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

void write_to_file(const std::vector<sliceable_set_t>& sets,
                   const std::filesystem::path& path) {
  constexpr std::size_t num_bytes =
      min_bytes_to_represent_bits(sliceable_set_t().size());
  std::ofstream file(path, std::ios::binary);
  for (const sliceable_set_t& set : sets) {
    char bytes[num_bytes];
    sliceable_set_to_bytes(set, bytes, num_bytes);
    file.write(bytes, num_bytes);
  }
}

char* read_from_file(const std::filesystem::path& path) {
  const std::size_t filesize = std::filesystem::file_size(path);
  char* array = new char[filesize];
  std::ifstream file(path, std::ios::binary);
  file.read(array, filesize);
  return array;
}

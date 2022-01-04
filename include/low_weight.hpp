#ifndef N_CUBE_LOW_WEIGHT_H_
#define N_CUBE_LOW_WEIGHT_H_

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "edge.hpp"
#include "prettyprint.hpp"
#include "sliceable_set.hpp"
#include "vertex.hpp"

/**
 *  Returns the slicing by a low weight halfspace.
 *
 *  The low weight halfspace is given by its normal vector and threshold
 *  (distance to the origin).
 **/
template <int32_t N>
sliceable_set_t<N> low_weight_halfspace_to_sliceable_set(
    const std::vector<int32_t>& normal, int32_t threshold,
    const std::vector<edge_t>& edges) {
  sliceable_set_t<N> ss;
  for (const auto& e : edges) {
    int32_t u_scalar = 0, v_scalar = 0;
    for (int32_t i = 0; i < N; ++i) {
      const auto u_i = get_coordinate(e.first, i);
      const auto v_i = get_coordinate(e.second, i);
      u_scalar += u_i * normal[i];
      v_scalar += v_i * normal[i];
    }
    if ((u_scalar < threshold && v_scalar > threshold) ||
        (u_scalar > threshold && v_scalar < threshold)) {
      ss[edge_to_int(e, edges)] = true;
    }
  }
  return ss;
}

/**
 *  Changes a normal vector containing only values in {-1, 1} into the next
 *  normal vector containing only values in {-1, 1}.
 *
 *  Returns true if the resulting normal vector is all -1 and false otherwise.
 *
 *  Naturally, the first call should be on a normal vector that is all -1.
 **/
bool next_one_weight_vector(std::vector<int32_t>& halfspace) {
  for (auto it = halfspace.rbegin(); it != halfspace.rend(); ++it) {
    if (*it == -1) {
      *it = 1;
      return true;
    } else {
      *it = -1;
    }
  }
  return false;
}

/**
 *  Changes a normal vector containing only values in {-max, ..., max} into the
 *  next normal vector containing only values in {-max, ..., max}.
 *
 *  Returns true if the resulting normal vector is all -max and false otherwise.
 *
 *  Naturally, the first call should be on a normal vector that is all -max.
 **/
bool next_low_weight_vector(std::vector<int32_t>& halfspace, int32_t max) {
  for (auto it = halfspace.rbegin(); it != halfspace.rend(); ++it) {
    if (*it == max) {
      *it = -max;
    } else {
      *it += 1;
      return true;
    }
  }
  return false;
}

/**
 *  Returns all slicings by low weight halfspaces satisfying the following:
 *    - The normal vector contains only values in {-1, 1}.
 *    - The threshold (distance to the origin) is one of the given thresholds.
 **/
template <int32_t N>
std::vector<sliceable_set_t<N>> compute_one_weight_sliceable_sets(
    const std::vector<int32_t>& thresholds, const std::vector<edge_t>& edges) {
  std::unordered_set<sliceable_set_t<N>> sets;
  std::vector<int32_t> normal(N, -1);
  do {
    for (const auto& threshold : thresholds) {
      const sliceable_set_t<N> ss =
          low_weight_halfspace_to_sliceable_set<N>(normal, threshold, edges);
      if (ss.any()) {
        sets.insert(ss);
      }
    }
  } while (next_one_weight_vector(normal));
  return std::vector<sliceable_set_t<N>>(sets.begin(), sets.end());
}

/**
 *  Returns all slicings by low weight halfspaces satisfying the following:
 *    - The normal vector contains only values in {-max, ..., max}.
 *    - The threshold (distance to the origin) is one of the given thresholds.
 **/
template <int32_t N>
std::vector<sliceable_set_t<N>> compute_low_weight_sliceable_sets(
    int32_t max, const std::vector<int32_t>& thresholds,
    const std::vector<edge_t>& edges) {
  std::unordered_set<sliceable_set_t<N>> sets;
  std::vector<int32_t> normal(N, -max);
  do {
    for (const auto& threshold : thresholds) {
      const sliceable_set_t<N> ss =
          low_weight_halfspace_to_sliceable_set<N>(normal, threshold, edges);
      if (ss.any()) {
        sets.insert(ss);
      }
    }
  } while (next_low_weight_vector(normal, max));
  return std::vector<sliceable_set_t<N>>(sets.begin(), sets.end());
}

/**
 *  Returns all maximal slicings by low weight halfspaces satisfying the
 *  following:
 *    - The normal vector contains only values in {-max, ..., max}.
 *    - The threshold (distance to the origin) is any integer value.
 **/
template <int32_t N>
std::vector<sliceable_set_t<N>> compute_low_weight_mss(
    const std::vector<edge_t>& edges, int32_t max) {
  std::vector<sliceable_set_t<N>> sets;
  std::vector<int32_t> normal(N, -max);
  do {
    for (int32_t d = 0; d < max * N; ++d) {
      const sliceable_set_t<N> ss =
          low_weight_halfspace_to_sliceable_set<N>(normal, d, edges);
      if (ss.any()) {
        const auto is_superset = [ss](const sliceable_set_t<N>& other) {
          return (other | ss) == other;
        };
        if (std::none_of(sets.begin(), sets.end(), is_superset)) {
          const auto is_subset = [ss](const sliceable_set_t<N>& other) {
            return (other | ss) == ss;
          };
          const auto it = remove_if(sets.begin(), sets.end(), is_subset);
          sets.erase(it, sets.end());
          sets.push_back(ss);
        }
      }
    }
  } while (next_low_weight_vector(normal, max));
  return sets;
}

/**
 *  Returns the smallest k so that the given sliceable sets contain a
 *  k-sliceable set.
 *
 *  This function uses significantly more memory than the below function which
 *  prevents efficient parallelization.
 **/
template <int32_t N>
int32_t combine_low_weight_sliceable_sets(
    const std::vector<sliceable_set_t<N>>& sets,
    const std::vector<edge_trans_t>& transformations) {
  std::vector<sliceable_set_t<N>> sets_k;
  for (const auto& ss : sets) {
    const auto usr = unique_sliceable_set<N>(ss, transformations);
    if (std::find(sets_k.begin(), sets_k.end(), usr) == sets_k.end()) {
      sets_k.push_back(usr);
    }
  }
  std::cout << "expanded size = " << sets.size() << std::endl;
  std::cout << "k = " << 1 << " size = " << sets_k.size() << std::endl;
  for (int i = 2;; ++i) {
    sets_k = combine_usr_mss<N>(sets_k, sets, transformations);
    std::cout << "k = " << i << " size = " << sets_k.size() << std::endl;
    for (const auto& ss : sets_k) {
      if (ss.all()) {
        return i;
      }
    }
  }
}

/**
 *  Returns the smallest k so that the given sliceable sets contain a
 *  k-sliceable set.
 **/
template <int32_t N>
int32_t combine_low_weight_sliceable_sets(
    const std::vector<sliceable_set_t<N>>& sets,
    const std::vector<edge_t>& edges) {
  std::vector<sliceable_set_t<N>> sets_k;
  for (const auto& ss : sets) {
    const auto usr = unique_sliceable_set<N>(ss, edges);
    if (std::find(sets_k.begin(), sets_k.end(), usr) == sets_k.end()) {
      sets_k.push_back(usr);
    }
  }
  std::cout << "expanded size = " << sets.size() << std::endl;
  std::cout << "k = " << 1 << " size = " << sets_k.size() << std::endl;
  for (int i = 2;; ++i) {
    sets_k = combine_usr_mss_parallel<N>(sets_k, sets, edges);
    std::cout << "k = " << i << " size = " << sets_k.size() << std::endl;
    for (const auto& ss : sets_k) {
      if (ss.all()) {
        return i;
      }
    }
  }
}

/**
 *  Writes in lexicographic order all slicings by low weight halfspaces
 *  satisfying the following to a file at the given path:
 *    - The normal vector contains only values in {-1, 1}.
 *    - The threshold (distance to the origin) is one of the given thresholds.
 *
 *  For each slicing the following is written in text:
 *    - The bitstring encoding of the slicing.
 *    - The normal vector.
 *    - The threshold.
 **/
template <int32_t N>
void write_one_weight_halfspaces_to_file(const std::vector<int32_t>& thresholds,
                                         const std::vector<edge_t>& edges,
                                         const std::filesystem::path& path) {
  std::vector<std::string> output;
  std::vector<int32_t> normal(N, -1);
  do {
    for (const auto& threshold : thresholds) {
      const sliceable_set_t<N> ss =
          low_weight_halfspace_to_sliceable_set<N>(normal, threshold, edges);
      if (ss.any()) {
        std::stringstream str_stream;
        str_stream << ss << " " << normal << " " << threshold;
        output.push_back(str_stream.str());
      }
    }
  } while (next_one_weight_vector(normal));
  std::sort(output.begin(), output.end());
  std::ofstream file(path);
  for (const auto& str : output) {
    file << str << std::endl;
  }
}

/**
 *  Writes in lexicographic order all slicings by low weight halfspaces
 *  satisfying the following to a file at the given path:
 *    - The normal vector contains only values in {-max, ..., max}.
 *    - The threshold (distance to the origin) is any integer value.
 *
 *  For each slicing the following is written in text:
 *    - The bitstring encoding of the slicing.
 *    - The normal vector.
 *    - The threshold.
 **/
template <int32_t N>
void write_low_weight_halfspaces_to_file(int32_t max,
                                         const std::vector<edge_t>& edges,
                                         const std::filesystem::path& path) {
  std::vector<std::string> output;
  std::vector<int32_t> normal(N, -max);
  do {
    for (int32_t threshold = 0; threshold < max * N; ++threshold) {
      const sliceable_set_t<N> ss =
          low_weight_halfspace_to_sliceable_set<N>(normal, threshold, edges);
      if (ss.any()) {
        std::stringstream str_stream;
        str_stream << ss << " " << normal << " " << threshold;
        output.push_back(str_stream.str());
      }
    }
  } while (next_low_weight_vector(normal, max));
  std::sort(output.begin(), output.end());
  std::ofstream file(path);
  for (const auto& str : output) {
    file << str << std::endl;
  }
}

#endif  // N_CUBE_LOW_WEIGHT_H_

#ifndef N_CUBE_MULTITHREADED_H_
#define N_CUBE_MULTITHREADED_H_

#include <algorithm>
#include <cstdint>
#include <thread>
#include <vector>

#include "bitset_comparator.hpp"
#include "edge.hpp"
#include "sliceable_set.hpp"
#include "vertex.hpp"

namespace ncube {

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
 *  Stores the unique symmetric representatives of the pairwise unions of two
 *  lists of sliceable sets in a range. Does NOT discard duplicates.
 **/
template <int32_t N>
void pairwise_unions_all(
    typename std::vector<sliceable_set_t<N>>::const_iterator sets_1_begin,
    typename std::vector<sliceable_set_t<N>>::const_iterator sets_1_end,
    typename std::vector<sliceable_set_t<N>>::const_iterator sets_2_begin,
    typename std::vector<sliceable_set_t<N>>::const_iterator sets_2_end,
    typename std::vector<sliceable_set_t<N>>::iterator unions_begin,
    const edge_lexicon_t<N>& edges) {
  auto unions_it = unions_begin;
  for (auto set_1 = sets_1_begin; set_1 != sets_1_end; ++set_1) {
    for (auto set_2 = sets_2_begin; set_2 != sets_2_end; ++set_2) {
      *unions_it = unique_sliceable_set<N>(*set_1 | *set_2, edges);
      ++unions_it;
    }
  }
}

/**
 *  Returns the unique symmetric representatives of the pairwise unions of two
 *  lists of sliceable sets.
 *
 *  This function is parallelized but requires a significant amount of memory.
 **/
template <int32_t N>
std::vector<sliceable_set_t<N>> pairwise_unions_parallel(
    const std::vector<sliceable_set_t<N>>& sets_1,
    const std::vector<sliceable_set_t<N>>& sets_2,
    const edge_lexicon_t<N>& edges) {
  // ensure effective parallelization
  if (sets_1.size() > sets_2.size()) {
    return pairwise_unions_parallel<N>(sets_2, sets_1, edges);
  }
  // divide mss
  const unsigned int num_threads = std::thread::hardware_concurrency();
  const auto set_1_workload = assign_workload(sets_2.size(), num_threads);
  // compute the unique symmetric representatives of all pairwise unions
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

}  // namespace ncube

#endif  // N_CUBE_MULTITHREADED_H_

#ifndef ALL_VS_ALL2_HPP
#define ALL_VS_ALL2_HPP
#include "utilities/matrix_partition.hpp"
#include <cstddef>
#include <vector>

template<class PairWiseComparator,class ContainerClass>
static inline void compare_pairs_in_range(bool triangle,
                                          PairWiseComparator *pw_comparator,
                                          const ContainerClass &references,
                                          size_t start_ref,size_t end_ref,
                                          const ContainerClass &queries,
                                          size_t start_query,size_t end_query)
{
#undef SHOW_COMPARISON_RANGES
#ifdef SHOW_COMPARISON_RANGES
  std::cerr << "# compare db=[" << start_ref << "," << end_ref
            << "] vs query=[" << start_query << "," << end_query
            << "]" << std::endl;
#endif
  for (size_t i = start_ref; i < end_ref; i++)
  {
    const size_t this_start_query = triangle ? i + 1 : start_query;
    pw_comparator->preprocess(i,references[i]);
    for (size_t j = this_start_query; j < end_query; j++)
    {
      if (pw_comparator->compare(i,j,queries[j]))
      {
        break;
      }
    }
  }
}

template<class PairWiseComparator,class ContainerClass>
static inline void all_against_all_compare_pairs
                     (size_t thread_id,
                      size_t task_num,
                      const ContainerClass &references,
                      const ContainerClass &queries,
                      bool same_container,
                      const MatrixPartition &matrix_partition,
                      const std::vector<PairWiseComparator *>
                        &pw_comparator_vector)
{
  const MatrixPartitionIntervalPair &itv = matrix_partition[task_num];
  const bool triangle = same_container && itv[2] == 0;

  compare_pairs_in_range<PairWiseComparator,ContainerClass>
                        (triangle,
                         pw_comparator_vector[thread_id],
                         references,
                         itv[0],
                         itv[0] + itv[1],
                         queries,
                         triangle ? itv[0] : itv[2],
                         triangle ? itv[0] + itv[1] : itv[2] + itv[3]);
}
#endif

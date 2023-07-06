#ifndef SORT_TRACK_PERMUTATION_HPP
#define SORT_TRACK_PERMUTATION_HPP
#include <cstddef>

template<typename T>
static inline bool sort_track_permutation(T *sortspace,
                                          size_t *permutation,
                                          const T *input,
                                          size_t num_elements)
{
  for (size_t idx = 0; idx < num_elements; idx++)
  {
    sortspace[idx] = input[idx];
    permutation[idx] = idx;
  }
  bool swapped = false;
  for (size_t pm = 0; pm < num_elements; pm++)
  {
    for (size_t pl = pm; pl > 0 and sortspace[pl-1] > sortspace[pl]; pl--)
    {
      const T tmp_cc = sortspace[pl-1];
      sortspace[pl-1] = sortspace[pl];
      sortspace[pl] = tmp_cc;
      const size_t tmp_idx = permutation[pl-1];
      permutation[pl-1] = permutation[pl];
      permutation[pl] = tmp_idx;
      swapped = true;
    }
  }
  return swapped;
}
#endif

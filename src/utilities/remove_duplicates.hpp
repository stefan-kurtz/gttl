#ifndef REMOVE_DUPLICATES_HPP
#define REMOVE_DUPLICATES_HPP
#include <cstddef>
#include <vector>

/* Function to remove duplicates in a sorted array. Does not use branches
   in inner loop. Returns number of pairwise distinct elements. */

template<typename T>
static size_t remove_duplicates(T *arr,size_t len)
{
  size_t n_dup = 0;
  for (size_t idx = 1; idx < len; idx++)
  {
    n_dup += (arr[idx-1-n_dup] == arr[idx]);
    arr[idx-n_dup] = arr[idx];
  }
  return len - n_dup;
}

/* When a vector is given we call this function with the data pointer
   and resize the vector afterwords. */

template<typename T>
static void remove_duplicates(std::vector<T> *vec)
{
  size_t n = remove_duplicates(vec->data(),vec->size());
  vec->resize(n);
}
#endif

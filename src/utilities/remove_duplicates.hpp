#ifndef REMOVE_DUPLICATES_HPP
#define REMOVE_DUPLICATES_HPP

/* a function to remove duplicates in a sorted array. Does not use branches
   in inner loop. Returns number of pairwise distinct elements. 
   is_duplicate must return 1 if elements are identical, otherwise 0 */

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
#endif

#ifndef BUCKETS_HPP
#define BUCKETS_HPP
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <utility>

/* A class to represent disjoint intervals of super interval beginning
   with index 0 such that all positions of the superinterval are covered
   by the intervals. These are specified by an array of their ends, where
   each end is excluded from the interval. The intervals are called buckets,
   hence the name. */

template<typename basetype>
class Buckets
{
  struct Iterator
  {
    private:
      size_t bucket_idx, num_buckets;
      basetype bucket_start;
      const basetype *bucket_ends;
    public:
      Iterator(const basetype *_bucket_ends,size_t _num_buckets) :
        bucket_idx(0),
        num_buckets(_num_buckets),
        bucket_start(0),
        bucket_ends(_bucket_ends) {}
     Iterator& operator++() /* prefix increment */
     {
       assert(bucket_idx < num_buckets);
       bucket_start = bucket_ends[bucket_idx++];
       return *this;
     }
     const std::pair<basetype,basetype> operator*(void) const
     {
       assert(bucket_idx < num_buckets);
       return {bucket_start,bucket_ends[bucket_idx]};
     }
     bool operator != (const Iterator& other) const noexcept
     {
       return bucket_idx != other.num_buckets;
     }
  };
  private:
    basetype *bucket_ends;
    size_t num_buckets;
  public:
  Buckets(size_t _num_buckets)
    : bucket_ends(new basetype [_num_buckets])
    , num_buckets(_num_buckets) {}
  ~Buckets(void)
  {
    delete[] bucket_ends;
  }
  void set(size_t idx,basetype value) noexcept
  {
    assert(idx < num_buckets);
    bucket_ends[idx] = value;
  }
  basetype *reference(void) noexcept
  {
    return bucket_ends;
  }
  [[nodiscard]] size_t size(void) const noexcept { return num_buckets; }
  [[nodiscard]] Iterator begin() const noexcept
  {
    return Iterator(bucket_ends,num_buckets);
  }
  [[nodiscard]] Iterator end() const { return Iterator(nullptr, num_buckets); }
  [[nodiscard]] basetype maximum_width_get(void) const noexcept
  {
    basetype bucket_start = 0;
    basetype maximum_bucket_width = 0;
    for (size_t idx = 0; idx < num_buckets; idx++)
    {
      assert(bucket_start <= bucket_ends[idx]);
      const basetype bucket_width = bucket_ends[idx] - bucket_start;
      maximum_bucket_width = std::max(maximum_bucket_width,bucket_width);
      bucket_start = bucket_ends[idx];
    }
    return maximum_bucket_width;
  }
};
#endif

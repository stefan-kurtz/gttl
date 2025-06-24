#ifndef MATRIX_PARTITION_HPP
#define MATRIX_PARTITION_HPP
#include <cstddef>
#include <cassert>
#include <cstdlib>
#include <array>
#include <vector>

using MatrixPartitionIntervalPair = std::array<size_t,4>;

static size_t matrix_partition_antidiagonal(const MatrixPartitionIntervalPair
                                              &mp)
{
  return mp[3] == 0 ? (mp[0] + mp[0]) : (mp[0] + mp[2]);
}

static size_t matrix_partition_antidiagonal(const MatrixPartitionIntervalPair
                                              *mp)
{
  return matrix_partition_antidiagonal(*mp);
}

class MatrixPartition
{
  using MatrixPartitionVector = std::vector<MatrixPartitionIntervalPair>;
  MatrixPartitionVector itv_list;
  MatrixPartitionIntervalPair split_interval(size_t a,size_t b) const noexcept
  {
    const size_t h = b/2 + (b % 2);
    return MatrixPartitionIntervalPair{a,h,a+h,b-h};
  }
  static int compare_itv(const void *va,const void *vb)
  {
    const MatrixPartitionIntervalPair *a
      = static_cast<const MatrixPartitionIntervalPair *>(va);
    const MatrixPartitionIntervalPair *b
      = static_cast<const MatrixPartitionIntervalPair *>(vb);
    const size_t anti_a = matrix_partition_antidiagonal(a);
    const size_t anti_b = matrix_partition_antidiagonal(b);
    if (anti_a < anti_b) { return -1; }
    if (anti_a > anti_b) { return 1; }
    if (std::get<0>(*a) < std::get<0>(*b)) { return -1; }
    if (std::get<0>(*a) > std::get<0>(*b)) { return +1; }
    assert(false);
    return 0;
  }

  public:
  MatrixPartition(size_t cutlen,size_t m,size_t n)
    : itv_list({})
  {
    assert(cutlen > 0);
    MatrixPartitionVector stack{};
    stack.push_back(MatrixPartitionIntervalPair{0,m,0,n});
    while (!stack.empty())
    {
      MatrixPartitionIntervalPair next_itv = stack.back();
      stack.pop_back();
      const size_t j = std::get<1>(next_itv);
      const size_t l = std::get<3>(next_itv);
      if (j <= cutlen && l <= cutlen)
      {
        itv_list.push_back(next_itv);
      } else
      {
        const size_t i = std::get<0>(next_itv);
        const size_t k = std::get<2>(next_itv);
        if (j < l)
        {
          MatrixPartitionIntervalPair itv = split_interval(k,l);
          stack.push_back({i,j,std::get<0>(itv),std::get<1>(itv)});
          stack.push_back({i,j,std::get<2>(itv),std::get<3>(itv)});
        } else
        {
          MatrixPartitionIntervalPair itv = split_interval(i,j);
          stack.push_back({std::get<0>(itv),std::get<1>(itv),k,l});
          stack.push_back({std::get<2>(itv),std::get<3>(itv),k,l});
        }
      }
    }
    qsort(itv_list.data(),itv_list.size(),sizeof(MatrixPartitionIntervalPair),
          compare_itv);
  }
  MatrixPartition(size_t cutlen,size_t m)
    : itv_list({})
  {
    assert(cutlen > 0);
    for (size_t idx = 0; idx < m; idx += cutlen)
    {
      if (idx + cutlen <= m)
      {
        itv_list.push_back(MatrixPartitionIntervalPair{idx,cutlen,0,0});
      } else
      {
        if (m > idx)
        {
          itv_list.push_back(MatrixPartitionIntervalPair{idx,m-idx,0,0});
        }
      }
    }
    const size_t num_simple_pairs = itv_list.size();
    assert(num_simple_pairs > 0);
    for (size_t i = 0; i < num_simple_pairs - 1; i++)
    {
      const MatrixPartitionIntervalPair itv_i = itv_list[i];
      for (size_t j = i+1; j < num_simple_pairs; j++)
      {
        const MatrixPartitionIntervalPair itv_j = itv_list[j];
        itv_list.push_back(MatrixPartitionIntervalPair{itv_i[0],itv_i[1],
                                                       itv_j[0],itv_j[1]});
      }
    }
    qsort(itv_list.data(),itv_list.size(),sizeof(MatrixPartitionIntervalPair),
          compare_itv);
  }
  MatrixPartitionIntervalPair operator [](size_t idx) const noexcept
  {
    assert(idx < itv_list.size());
    return itv_list[idx];
  }
  using ConstIterator = MatrixPartitionVector::const_iterator;
  ConstIterator begin(void) const
  {
    return itv_list.cbegin();
  }

  ConstIterator end(void) const
  {
    return itv_list.cend();
  }
  size_t size(void) const noexcept
  {
    return itv_list.size();
  }
};
#endif

#ifndef MATRIX_PARTITION_HPP
#define MATRIX_PARTITION_HPP
#include <cstddef>
#include <algorithm>
#include <array>

using MatrixPartitionIntervalPair = std::array<size_t,4>;

class MatrixPartition
{
  private:
  MatrixPartitionIntervalPair split_interval(size_t a,size_t b)
  {
    const size_t h = b/2 + (b % 2);
    return MatrixPartitionIntervalPair{a,h,a+h,b-h};
  }
  std::vector<MatrixPartitionIntervalPair> itv_list{};

  public:
  MatrixPartition(size_t cutlen,size_t m,size_t n)
  {
    std::vector<MatrixPartitionIntervalPair> stack{};
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
  }
  const std::vector<MatrixPartitionIntervalPair> &intervals(void) const noexcept
  {
    return itv_list;
  }
};
#endif

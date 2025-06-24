#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include "utilities/matrix_partition.hpp"

int main(int argc,char *argv[])
{
  int cutlen;
  int rows;
  int cols;
  if (argc != 4 ||
      sscanf(argv[1],"%d",&cutlen) != 1 || cutlen < 1 ||
      sscanf(argv[2],"%d",&rows) != 1 || rows < 1 ||
      sscanf(argv[3],"%d",&cols) != 1 || cols < 0)
  {
    std::cerr << "Usage: " << argv[0] << " <cutlen> <rows> <cols>\n";
    return EXIT_FAILURE;
  }
  MatrixPartition mp = cols == 0 ? MatrixPartition(static_cast<size_t>(cutlen),
                                                   static_cast<size_t>(rows))
                                 : MatrixPartition(static_cast<size_t>(cutlen),
                                                   static_cast<size_t>(rows),
                                                   static_cast<size_t>(cols));
  size_t num_pairs = 0;
  for (auto && itv : mp)
  {
    std::cout << "((" << itv[0] << ", " << itv[1] << "), ";
    size_t this_pairs;
    if (itv[3] == 0)
    {
      this_pairs = (itv[1] * (itv[1] - 1))/2;
      std::cout << "None" << ")\t" << matrix_partition_antidiagonal(itv) << "\t"
                << this_pairs << '\n';
    } else
    {
      this_pairs = itv[1] * itv[3];
      std::cout << "(" << itv[2] << ", " << itv[3] << "))\t"
                << matrix_partition_antidiagonal(itv) << "\t" << this_pairs
                << '\n';
    }
    num_pairs += this_pairs;
  }
  std::cout << "# number of parts\t" << mp.size() << '\n';
#ifndef NDEBUG
  size_t expected;
#endif
  if (cols == 0)
  {
#ifndef NDEBUG
    const size_t n = static_cast<size_t>(rows);
    expected = (n * (n-1))/2;
#endif
  } else
  {
#ifndef NDEBUG
    expected = static_cast<size_t>(rows) * static_cast<size_t>(cols);
#endif
  }
  assert(num_pairs == expected);
  std::cout << "# number of pairs\t" << num_pairs << '\n';
  return EXIT_SUCCESS;
}

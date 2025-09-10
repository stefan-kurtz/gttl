#include <bitset>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "utilities/multibitvector.hpp"

int main(int argc, char* argv[])
{
  if (argc != 1)
  {
    std::cerr << "USAGE: ./" << argv[0] << '\n';
    return EXIT_FAILURE;
  }
  const std::vector<size_t> sizes = {1, 7, 63, 64, 65, 127, 128, 129};

  for (const size_t n_bits : sizes)
  {
    std::cout << "Testing n_bits = " << n_bits << '\n';

    Multibitvector<false> bv(n_bits);

    for (size_t i = 0; i < n_bits; ++i)
    {
      bv.set(i);
      if (not bv[i])
      {
        std::cerr << argv[0] << ": set(" << i << ") did not set bit!\n";
	return EXIT_FAILURE;
      }
      bv.reset(i);
      if (bv[i])
      {
        std::cerr << argv[0] << ": reset(" << i << ") did not reset bit!\n";
	return EXIT_FAILURE;
      }
    }

    Multibitvector<false> bv2(n_bits);
    size_t manual_count = 0;

    for (size_t i = 0; i < n_bits; i += 2)
    {
      bv2.set(i);
      ++manual_count;
    }
    const size_t c = bv2.count();
    if (c != manual_count)
    {
      std::cerr << argv[0] << ": count() mismatch. Expected " << manual_count << " got " << c << '\n';
      return EXIT_FAILURE;
    }

    Multibitvector<false> bv3(n_bits);
    bv3 |= bv2;
    if (bv3.count() != bv2.count())
    {
      std::cerr << argv[0] << ": |= mismatch from initially empty bitvector!\n";
      return EXIT_FAILURE;
    }

    std::bitset<256> golden;
    Multibitvector<false> bv4(n_bits);
    for (size_t i = 0; i < n_bits; ++i)
    {
      if (i % 3 == 0)
      {
        bv4.set(i);
	golden.set(i);
      }
    }
    for (size_t i = 0; i < n_bits; ++i)
    {
      if (bv4[i] != golden.test(i))
      {
        std::cerr << argv[0] << ": mismatch against std::bitset at bit " << i << '\n';
	return EXIT_FAILURE;
      }
    }
    if (bv4.count() != static_cast<size_t>(golden.count()))
    {
      std::cerr << argv[0] << ": count() mismatch against std::bitset!\n";
      return EXIT_FAILURE;
    }

    size_t expected = 0;
    Multibitvector<false> bv5(n_bits);
    for (size_t i = 0; i < n_bits; ++i)
    {
      bv5.set(i);
      ++expected;
    }
    if (bv5.count() != expected)
    {
      std::cerr << argv[0] << ": stray bits counted for size " << n_bits << '\n';
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}

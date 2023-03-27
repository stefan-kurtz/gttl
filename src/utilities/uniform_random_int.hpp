#ifndef UNIFORM_RANDOM_INT_HPP
#define UNIFORM_RANDOM_INT_HPP
#include <random>

template<typename IntType>
class UniformRandomInteger
{
  std::random_device seed;
  std::mt19937 generator;
  std::uniform_int_distribution<IntType> distribution;

  public:
  UniformRandomInteger(IntType low, IntType high, unsigned int own_seed = 0)
    : generator(own_seed == 0 ? seed() : own_seed)
    , distribution(low,high)
  {}
  IntType get(void)
  {
    return distribution(generator);
  }
};
#endif

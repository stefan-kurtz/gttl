#ifndef UNIFORM_RANDOM_INT_HPP
#define UNIFORM_RANDOM_INT_HPP
#include <cfloat>
#include <random>

template<typename IntType>
class UniformRandomInteger
{
  std::random_device seed{};
  std::mt19937 generator;
  std::uniform_int_distribution<IntType> distribution;

  void set(IntType low, IntType high)
  {
    std::uniform_int_distribution<IntType>::param_type param(low, high);
    distribution.param(param);
  }
  public:
  UniformRandomInt(IntType low, IntType high, unsigned int own_seed = 0)
    : generator(own_seed == 0 ? seed() : own_seed)
  {
    set(low, high);
  }
  IntType get(void)
  {
    return distribution(generator);
  }
};
#endif

#ifndef UNIFORM_RANDOM_DOUBLE_HPP
#define UNIFORM_RANDOM_DOUBLE_HPP
#include <cfloat>
#include <random>

class UniformRandomDouble
{
  std::random_device seed_gen; // generate a seed
  std::mt19937_64 generator; // mersenne_twister_engine
  std::uniform_real_distribution<double> distribution;

  public:
  UniformRandomDouble(double low = DBL_MIN,
                      double high = DBL_MAX,
                      unsigned int own_seed = 0)
    : generator(own_seed == 0 ? seed_gen() : own_seed)
    , distribution(low,high)
  {}
  double get(void)
  {
    return distribution(generator);
  }
};
#endif

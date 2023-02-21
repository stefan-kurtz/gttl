#ifndef UNIFORM_RANDOM_DOUBLE_HPP
#define UNIFORM_RANDOM_DOUBLE_HPP
#include <cfloat>
#include <random>

class UniformRandomDouble
{
  std::random_device seed{};
  std::mt19937 generator;
  std::uniform_real_distribution<double> distribution;

  public:
  UniformRandomDouble(double low = DBL_MIN,
                      double high = DBL_MAX,
                      unsigned int own_seed = 0)
    : generator(own_seed == 0 ? seed() : own_seed)
  {
    std::uniform_real_distribution<double>::param_type param(low, high);
    distribution.param(param);
  }
  double get(void)
  {
    return distribution(generator);
  }
};
#endif

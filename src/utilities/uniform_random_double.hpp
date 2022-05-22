#ifndef UNIFORM_RANDOM_DOUBLE_HPP
#define UNIFORM_RANDOM_DOUBLE_HPP
#include <cfloat>
#include <random>

class UniformRandomDouble
{
  std::random_device seed{};
  std::mt19937 generator;
  std::uniform_real_distribution<double> distribution;

  void set(double low, double high)
  {
    std::uniform_real_distribution<double>::param_type param(low, high);
    distribution.param(param);
  }
  public:
  UniformRandomDouble(unsigned int own_seed = 0)
    : generator(own_seed == 0 ? seed() : own_seed)
  {
    set(DBL_MIN, DBL_MAX);
  }
  UniformRandomDouble(double low, double high, unsigned int own_seed = 0)
    : generator(own_seed == 0 ? seed() : own_seed)
  {
    set(low, high);
  }
  double get(void)
  {
    return distribution(generator);
  }
};
#endif

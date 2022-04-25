#ifndef UNIFORM_RANDOM_DOUBLE_HPP
#define UNIFORM_RANDOM_DOUBLE_HPP
#include <cfloat>
#include <random>

class UniformRandomDouble
{
  std::random_device _rd{};
  std::mt19937 _gen{_rd()};
  std::uniform_real_distribution<double> _dist;

  void set(double low, double high)
  {
    std::uniform_real_distribution<double>::param_type param(low, high);
    _dist.param(param);
  }
  public:
  UniformRandomDouble(void)
  {
    set(DBL_MIN, DBL_MAX);
  }
  UniformRandomDouble(double low, double high)
  {
    set(low, high);
  }
  double get(void)
  {
    return _dist(_gen);
  }
};
#endif

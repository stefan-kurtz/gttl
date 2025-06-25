#ifndef NTTABLE_HPP
#define NTTABLE_HPP

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>
#include <limits>

class NtTableResult
{
  private:
  double const *f_n_ptr;
  const size_t t_max;
  const double F0;
  const double F1;
  const size_t F1_count;

  public:
  NtTableResult(double const *_f_n,
                size_t _t_max,
                double _F0,
                double _F1,
                size_t _F1_count)
    : f_n_ptr(_f_n)
    , t_max(_t_max)
    , F0(_F0)
    , F1(_F1)
    , F1_count(_F1_count)
  { }
  ~NtTableResult(void)
  {
    delete[] f_n_ptr;
  }
  double f_n_at(size_t idx) const
  {
    assert(idx <= t_max);
    return f_n_ptr[idx];
  }
  double F0_get(void) const { return F0; }
  double F1_get(void) const { return F1; }
  size_t F1_count_get(void) const { return F1_count; }
  size_t t_max_get(void) const { return t_max; }
  double count_estimate_get(void) const
  {
    return F0_get() + (F0_get() - f_n_at(1)) * 2.0;
  }
};

class NtTable
{
  private:
  using CountType = uint16_t;
  static constexpr const size_t CountTypeMax
    = static_cast<size_t>(std::numeric_limits<CountType>::max());
  const size_t s_value;
  const size_t r_value;
  const uint64_t s_mask;
  const uint64_t r_mask;
  size_t qgram_count,
         sequences_number;
  std::vector<CountType> table;

 public:
  NtTable(size_t _s, size_t _r)
    : s_value(_s)
    , r_value(_r)
    , s_mask(_s == 0 ? 0 : ((~uint64_t(0)) << (64 - _s)))
    , r_mask(               (~uint64_t(0)) >> (64 - _r))
    , qgram_count(0)
    , sequences_number(0)
    , table(size_t(1) << _r, 0)
  { }

  void add_hash(uint64_t hash)
  {
    qgram_count++;
    if ((hash & s_mask) != 0)
    {
      return;
    }
    table[hash & r_mask] += (table[hash & r_mask] < CountTypeMax);
  }

  void merge(const NtTable &other)
  {
    assert(s_value == other.s_value and r_value == other.r_value and
           table.size() == other.table.size());

    for (size_t idx = 0; idx < table.size(); idx++)
    {
      if (table[idx] < CountTypeMax - other.table[idx])
      {
        table[idx] += other.table[idx];
      } else // overflow
      {
        table[idx] = CountTypeMax;
      }
    }
    qgram_count += other.qgram_count;
    sequences_number += other.sequences_number;
  }
  double estimate_F0(void) const
  {
    size_t p0 = 0;
    for (auto t : table)
    {
      if (t == 0)
      {
        p0++;
      }
    }

    if (p0 == 0)
    {
      throw std::domain_error("The count is zero resulting in log(0).");
    }
    return -std::log(static_cast<double>(p0) / (uint64_t(1) << r_value)) *
           (uint64_t(1) << (s_value + r_value));
  }

  NtTableResult estimate_all(void) const
  {
    size_t t_max = 0;
    for (auto t : table)
    {
      if (t > t_max)
      {
        t_max = t;
      }
    }

    std::vector<size_t> p(t_max + 1, 0);
    for (auto t : table)
    {
      assert(t <= t_max);
      p[t]++;
    }

    if (p[0] == 0)
    {
      throw std::domain_error("a count of zero leads to undefined values");
    }

    static_assert(sizeof(double) == sizeof(size_t));
    // reuse memory for p for p_n, as the latter is just used for scaling
    double *const p_n = reinterpret_cast<double *>(p.data());
    for (size_t i = 0; i < p.size(); i++)
    {
      p_n[i] = static_cast<double>(p[i]) / (size_t(1) << r_value);
    }

    const double F0 = -std::log(p_n[0]) * (uint64_t(1) << (s_value + r_value));
    // allocate memory here and pass pointer to NtTableResult where it
    // is deleted once NtTableResult gets out of scope
    double *const f_n = new double[t_max + 1];
    f_n[0] = 0;
    const double denominator = static_cast<double>(1.0) /
                               (p_n[0] * std::log(p_n[0]));
    for (size_t i = 1; i <= t_max; i++)
    {
      double sum = 0.0;
      for (size_t j = 1; j < i; j++)
      {
        sum += static_cast<double>(j) * p_n[i - j] * f_n[j];
      }

      f_n[i] = -p_n[i] * denominator -
               (1.0 / (static_cast<double>(i) * p_n[0])) * sum;
    }

    double F1 = 0.0;
    for (size_t i = 1; i <= t_max; i++)
    {
      f_n[i] *= F0;
      F1 += static_cast<double>(i) * f_n[i];
    }

    return NtTableResult(f_n, t_max, F0, F1, qgram_count);
  }

  size_t F1_count_get(void) const noexcept
  {
    return qgram_count;
  }

  void sequences_number_set(size_t _sequences_number)
  {
    assert(sequences_number == 0);
    sequences_number = _sequences_number;
  }

  size_t sequences_number_get(void) const noexcept
  {
    return sequences_number;
  }
};
#endif

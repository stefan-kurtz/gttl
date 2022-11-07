#ifndef BLAST_STAT_HPP
#define BLAST_STAT_HPP
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cassert>
#include <cstdbool>
#include "utilities/str_format.hpp"

class BlastStatistics
{
  using BlastStatisticsLine = struct
  {
    int gap_open_penalty,
        gap_extension_penalty;
    double lambda, kappa;
  };
  static constexpr const BlastStatisticsLine blosum62_stat[]
  {
    {11, 2, 0.297,  0.082},
    {10, 2, 0.291,  0.075},
    { 9, 2, 0.279,  0.058},
    { 8, 2, 0.264,  0.045},
    { 7, 2, 0.239,  0.027},
    { 6, 2, 0.201,  0.012},
    {13, 1, 0.292,  0.071},
    {12, 1, 0.283,  0.059},
    {11, 1, 0.267,  0.041},
    {10, 1, 0.243,  0.024},
    { 9, 1, 0.206,  0.010}
  };
  static constexpr const BlastStatisticsLine blosum62_scaled_stat[] =
  {
    {44, 4, 0.08354, 0.08526}
  };
  int blast_stat_available(const BlastStatisticsLine *stat,
                           const size_t items,
                           int8_t gap_open_penalty,
                           int8_t gap_extension_penalty)
  {
    for (size_t idx = 0; idx < items; idx++)
    {
      if (stat[idx].gap_open_penalty == gap_open_penalty &&
          stat[idx].gap_extension_penalty == gap_extension_penalty)
      {
        return idx;
      }
    }
    return -1;
  }
  double log_kappa_d_log2, lambda_d_log2;

  public:
  BlastStatistics (int8_t gap_open_penalty,
                   int8_t gap_extension_penalty,
                   bool scaled)
  {
    const BlastStatisticsLine *stat = scaled ? &blosum62_scaled_stat[0]
                                             : &blosum62_stat[0];
    const size_t items
      = scaled ? sizeof blosum62_scaled_stat/sizeof blosum62_scaled_stat[0]
               : sizeof blosum62_stat/sizeof blosum62_stat[0];
    const int idx = blast_stat_available(stat,items,gap_open_penalty,
                                         gap_extension_penalty);
    if (idx == -1)
    {
      StrFormat msg(": no Gumbel parameters for computing bits scores "
                    "available for blosum62 matrix and gap parameters %d/%d",
                    gap_open_penalty,
                    gap_extension_penalty);
      throw msg.str();
    }
    assert(idx >= 0 && idx < static_cast<int>(items));
    const double lambda = stat[idx].lambda,
                 kappa = stat[idx].kappa;
    log_kappa_d_log2 = log(kappa)/log(2.0);
    lambda_d_log2 = lambda/log(2.0);
#ifdef GUMBLE_OUT
    printf("lambda_d_log2=%.6e\n",blast_stat->lambda_d_log2);
#endif
  }
  uint32_t raw_score2bit_score(uint32_t raw_score) const noexcept
  {
    return static_cast<uint32_t>(floor(lambda_d_log2 *
                                       static_cast<double>(raw_score) -
                                       log_kappa_d_log2 + 0.5));
  }
};
#endif

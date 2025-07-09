#ifndef ALIGN_POLISH_HPP
#define ALIGN_POLISH_HPP
#include <algorithm>
#include <cstddef>
#include <cassert>
#include <string>
#include <cstdint>
#include <iostream>

class AlignmentPolishing
{
  using AlignmentPolishingValue = struct { int16_t score_sum, diff_from_max;};
  private:
  size_t lag_last_d_with_pp,
         cut_depth,
         num_entries,
         history_mask;
  int64_t match_score,
          difference_score;
  AlignmentPolishingValue *values;
  void fill_polishing_info(size_t currentdepth,
                           size_t prefix,
                           int64_t score,
                           int64_t maxscore)
  {
    assert(currentdepth <= cut_depth);
    if (currentdepth == cut_depth)
    {
      assert(prefix < num_entries && score >= INT16_MIN + maxscore);
      values[prefix].diff_from_max = static_cast<int16_t>(score - maxscore);
      values[prefix].score_sum = static_cast<int16_t>(score);
    } else
    {
      if (score > maxscore)
      {
        maxscore = score;
      }
      assert(score >= INT16_MIN - difference_score);
      fill_polishing_info(currentdepth+1, prefix << 1,
                          score - difference_score,maxscore);
      assert(score <= INT16_MAX - match_score);
      fill_polishing_info(currentdepth+1,(prefix << 1) | size_t(1),
                          score + match_score,maxscore);
    }
  }
  [[nodiscard]] std::string
  polish_intbits2string(size_t bits, size_t bs) const noexcept
  {
    std::string cs{};

    assert(bits > 0 && bits <= sizeof (size_t) * CHAR_BIT);
    for (size_t mask = (size_t(1)) << (bits-1); mask > 0; mask >>= 1)
    {
       cs.push_back((bs & mask) ? '1' : '0');
    }
    return cs;
  }
  size_t lag_last_d_with_pp_evaluate(double _polishing_error_percentage)
  {
    if (_polishing_error_percentage <= 2.0)
    {
      return 1;
    }
    if (_polishing_error_percentage <= 4.0)
    {
      return 2;
    }
    if (_polishing_error_percentage <= 8.0)
    {
      return 4;
    }
    return 5;
  }
  public:
  AlignmentPolishing(double _polishing_error_percentage,
                     size_t _history_size,
                     double _matchscore_bias)
    : lag_last_d_with_pp(lag_last_d_with_pp_evaluate(
                            _polishing_error_percentage))
    , cut_depth(_history_size == 0 ? size_t(15)
                                   : std::min(_history_size/2,size_t(15)))
    , num_entries(size_t(1) << cut_depth)
    , history_mask(num_entries - 1)
    , match_score(20.0 * _polishing_error_percentage * _matchscore_bias)
    , difference_score(1000.0 - match_score)
    , values(new AlignmentPolishingValue [num_entries])
  {
    assert(match_score <= 1000.0);
    fill_polishing_info(0,0,0,0);
  }
  ~AlignmentPolishing(void)
  {
    delete[] values;
  }

  [[nodiscard]] bool
  is_polished_brute_force(uint64_t matchhistory, bool withoutput) const noexcept
  {
    size_t idx;
    uint64_t mask;
    int64_t sum_score = 0;

    for (mask = uint64_t(1), idx = 0; idx < 2 * cut_depth; idx++, mask <<= 1)
    {
      if (matchhistory & mask)
      {
        sum_score += match_score;
      } else
      {
        sum_score -= difference_score;
      }
      if (withoutput)
      {
        std::cout << idx << "\t" << sum_score << '\n';
      }
      if (sum_score < 0)
      {
        return false;
      }
    }
    return true;
  }
  [[nodiscard]] bool is_polished(uint64_t match_history) const noexcept
  {
    auto value0 = values[match_history & history_mask];
    auto value1 = values[(match_history >> cut_depth) & history_mask];
    return value0.diff_from_max >= 0 &&
           (value0.score_sum + value1.diff_from_max) >= 0;
  }
  [[nodiscard]] size_t lag_last_d_with_pp_get(void) const noexcept
  {
    return lag_last_d_with_pp;
  }
  void show(void) const noexcept
  {
    std::cout << "# cut_depth\t" << cut_depth << '\n';
    std::cout << "# entries\t" << num_entries << '\n';
    std::cout << "# match_score\t" << match_score << '\n';
    std::cout << "# difference_score\t" << difference_score << '\n';
    std::cout << "# lag_last_d_with_pp\t" << lag_last_d_with_pp << '\n';
    std::cout << "# history_mask\t"
              << polish_intbits2string(cut_depth, history_mask) << '\n';
    for (size_t idx = 0; idx < num_entries; idx++)
    {
      std::cout << "# " << polish_intbits2string(cut_depth, idx) << "\t"
                << values[idx].score_sum << "\t" << values[idx].diff_from_max
                << "\t"
                << (is_polished(static_cast<uint64_t>(idx)) ? "true" : "false")
                << '\n';
    }
  }
};
#endif

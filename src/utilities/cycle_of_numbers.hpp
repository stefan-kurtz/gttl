#ifndef CYCLE_OF_NUMBERS_HPP
#define CYCLE_OF_NUMBERS_HPP
#include <cstddef>
#include <cstdint>
#include <cassert>
#include <vector>

class CycleOfNumbers
{
  uint8_t current;
  std::vector<uint8_t> next_tab;

  public:

  CycleOfNumbers(const std::vector<uint8_t> &forbidden_characters)
    : current(0)
    , next_tab(UINT8_MAX+1,0)
  {
    for (size_t idx = 0; idx < forbidden_characters.size(); idx++)
    {
      const uint8_t this_forbidden_character = forbidden_characters[idx];
      /* character must be supplied in alphabetical order */
      assert(idx == 0 || forbidden_characters[idx-1] <
                         this_forbidden_character);
      next_tab[static_cast<int>(this_forbidden_character)] = uint8_t(1);
    }
    int next_value = -1;
    for (int idx = 0; idx <= UINT8_MAX; idx++)
    {
      if (next_tab[idx] == 0)
      {
        next_value = idx;
        break;
      }
    }
    assert(next_value != -1);
    for (int idx = static_cast<int>(UINT8_MAX); idx >= 0; idx--)
    {
      if (next_tab[idx] == 0)
      {
        assert(next_value >= 0 && next_value <= UINT8_MAX);
        next_tab[idx] = static_cast<uint8_t>(next_value);
        next_value = idx;
      }
    }
  }
  uint8_t next(void)
  {
    const uint8_t ret_value = current;
    current = next_tab[current];
    return ret_value;
  }
};
#endif

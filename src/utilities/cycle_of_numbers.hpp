#ifndef CYCLE_OF_NUMBERS_HPP
#define CYCLE_OF_NUMBERS_HPP
#include <cstdint>
#include <cassert>
#include <vector>

class CycleOfNumbers
{
  uint8_t current;
  std::vector<uint8_t> next_tab;

  public:

  CycleOfNumbers(const char *forbidden)
    : current(0)
    , next_tab(UINT8_MAX+1,0)
  {
    for (const char *f = forbidden; *f != '\0'; f++)
    {
      next_tab[static_cast<int>(*f)] = uint8_t(1);
    }
    assert(next_tab[UINT8_MAX] == 0 && next_tab[0] == 0);
    uint8_t next_value = UINT8_MAX;
    for (uint8_t idx = UINT8_MAX; idx > 0; idx--)
    {
      if (next_tab[idx-1] == 0)
      {
        assert(idx-1 < next_value);
        next_tab[idx-1] = next_value;
        next_value = idx-1;
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

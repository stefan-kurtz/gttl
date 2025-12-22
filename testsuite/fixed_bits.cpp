#include <cstdint>
#include <string>
#include <climits>
#include <iostream>
#include <cstring>
#include <format>
#include <cassert>

class FixedBitClass
{
  public:
  static constexpr const int num_bits = 36;
  private:
  static constexpr const int num_bytes = (2 * num_bits + CHAR_BIT - 1)/CHAR_BIT;
  struct TwoValues
  {
    uint64_t a:num_bits,
             b:num_bits;
  };
  union
  {
    uint8_t vector[num_bytes];
    TwoValues two_values;
  } overlay;
  public:
  FixedBitClass(uint64_t _a, uint64_t _b)
  {
#ifndef NDEBUG
    constexpr const uint64_t max_value = (uint64_t(1) << num_bits) - 1;
#endif
    assert(_a <= max_value and _b <= max_value);
    overlay.two_values.a = _a;
    overlay.two_values.b = _b;
  }
  [[nodiscard]] std::string to_string(void) const
  {
    std::string s = "(" + std::to_string(overlay.two_values.a) +
                    "," + std::to_string(overlay.two_values.b) +
                    ")";
    for (unsigned char idx : overlay.vector)
    {
      s += std::format(" {:b}",idx);
    }
    return s;
  }
  bool operator <(const FixedBitClass &other) const
  {
    return memcmp(&overlay.vector[0],&other.overlay.vector[0],num_bytes);
  }
};

int main(void)
{
  const uint64_t max_value = (uint64_t(1) << FixedBitClass::num_bits);
  const FixedBitClass fixed_bits(max_value, max_value);
  std::cout << fixed_bits.to_string() << '\n';
};

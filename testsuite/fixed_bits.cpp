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
  static constexpr const int num_bits0 = 32;
  static constexpr const int num_bits1 = 40;
  static constexpr const int num_bytes = (num_bits0 + num_bits1 + CHAR_BIT - 1)
                                         /CHAR_BIT;
  static constexpr const uint64_t max_value0 = (uint64_t(1) << num_bits0) - 1;
  static constexpr const uint64_t max_value1 = (uint64_t(1) << num_bits1) - 1;
  private:
  struct TwoValues
  {
    uint64_t a:num_bits0;
    uint64_t b:num_bits1;
  };
  union
  {
    uint8_t vector[num_bytes];
    TwoValues two_values;
  } overlay;
  public:
  FixedBitClass(uint64_t _a, uint64_t _b)
  {
    assert(_a <= max_value0 and _b <= max_value1);
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
      s += std::format(" {:03d}",idx);
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
  std::cout << "num_bytes=" << FixedBitClass::num_bytes << '\n';
  const FixedBitClass fixed_bits(FixedBitClass::max_value0,
                                 FixedBitClass::max_value1);
  std::cout << fixed_bits.to_string() << '\n';
};

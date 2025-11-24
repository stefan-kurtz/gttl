#ifndef BYTES_UNIT_HPP
#define BYTES_UNIT_HPP
#include <array>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <type_traits>
#include "utilities/constexpr_for.hpp"
#include "utilities/bitpacker.hpp"

template<size_t sizeof_unit,int bit_groups>
class BytesUnit
{
  using basetype = std::conditional_t<
                               sizeof_unit >= 8,
                               uint64_t,
                               std::conditional_t<sizeof_unit >= 4,
                                                  uint32_t,
                                                  uint16_t>>;

  private:
    uint8_t bytes[sizeof_unit];
  public:
    BytesUnit() = default;
    BytesUnit(const GttlBitPacker<sizeof_unit,bit_groups> &bitpacker,
              const std::array<uint64_t, bit_groups> &to_be_encoded)
    {
      static_assert(sizeof *this == sizeof_unit);
      static constexpr const int last_idx
        = sizeof_unit == sizeof(basetype) ? (bit_groups-1) : (bit_groups-2);
      assert(to_be_encoded[0] <= bitpacker.mask_tab[0]);
      basetype integer = to_be_encoded[0] << bitpacker.shift_tab[0];
      for (int idx = 1; idx <= last_idx; idx++)
      {
        assert(to_be_encoded[idx] <= bitpacker.mask_tab[idx]);
        integer |= (to_be_encoded[idx] << bitpacker.shift_tab[idx]);
      }
      if constexpr (sizeof_unit > sizeof(basetype))
      {
        const basetype overflow_value = to_be_encoded[bit_groups-1];
        assert(overflow_value <= bitpacker.max_overflow);
        integer |= (overflow_value >> bitpacker.overflow_left_shift);
        constexpr const int overflow = sizeof_unit -
                                       static_cast<int>(sizeof(basetype));
        constexpr_for<0,overflow,1>([&](auto this_idx)
        {
          constexpr const int shift = 8 * (overflow - this_idx - 1);
          bytes[sizeof(basetype)+this_idx]
            = static_cast<uint8_t>(overflow_value >> shift);
        });
#ifdef WITH_CASES
        if constexpr (sizeof_unit == sizeof(basetype) + 1)
        {
          bytes[sizeof(basetype)] = static_cast<uint8_t>(overflow_value);
        } else
        {
          if constexpr (sizeof_unit == sizeof(basetype) + 2)
          {
            bytes[sizeof(basetype)] = static_cast<uint8_t>(overflow_value >> 8);
            bytes[sizeof(basetype)+1] = static_cast<uint8_t>(overflow_value);
          } else
          {
            if constexpr (sizeof_unit == sizeof(basetype) + 3)
            {
              bytes[sizeof(basetype)]
                = static_cast<uint8_t>(overflow_value >> 16);
              bytes[sizeof(basetype)+1]
                = static_cast<uint8_t>(overflow_value >> 8);
              bytes[sizeof(basetype)+2]
                = static_cast<uint8_t>(overflow_value);
            } else
            {
              static_assert(sizeof_unit == sizeof(basetype) + 4);
              bytes[sizeof(basetype)]
                = static_cast<uint8_t>(overflow_value >> 24);
              bytes[sizeof(basetype)+1]
                = static_cast<uint8_t>(overflow_value >> 16);
              bytes[sizeof(basetype)+2]
                = static_cast<uint8_t>(overflow_value >> 8);
              bytes[sizeof(basetype)+3]
                = static_cast<uint8_t>(overflow_value);
            }
          }
        }
#endif
      }
      memcpy(&bytes[0],reinterpret_cast<uint8_t *>(&integer),sizeof integer);
    }

    template <int idx>
    [[nodiscard]] uint64_t
    decode_at(const GttlBitPacker<sizeof_unit, bit_groups> &bitpacker)
                                 const noexcept
    {
      static_assert(idx >= 0 && idx < bit_groups);
      static_assert(std::is_trivially_copyable_v<basetype>);
      basetype integer;
      std::memcpy(&integer, bytes, sizeof(basetype));

      if constexpr (sizeof_unit == sizeof(basetype) || idx < bit_groups - 1)
      {
        return static_cast<uint64_t>(integer >> bitpacker.shift_tab[idx]) &
               bitpacker.mask_tab[idx];
      } else /* sizeof_unit > sizeof(basetype) && idx == bit_groups - 1 */
      {
        constexpr const int overflow = sizeof_unit -
                                       static_cast<int>(sizeof(basetype));
        uint64_t this_value = integer << bitpacker.overflow_left_shift;
        constexpr_for<0,overflow,1>([&](auto this_idx)
        {
          constexpr const int shift = 8 * (overflow - this_idx - 1);
          this_value |= (static_cast<uint64_t>(bytes[sizeof(basetype)+this_idx])
                                               << shift);
        });
        return this_value & bitpacker.max_overflow;
#ifdef WITH_CASES
        if constexpr (sizeof_unit == sizeof(basetype) + 1)
        {
          return ((integer << bitpacker.overflow_left_shift) |
                  static_cast<uint64_t>(bytes[sizeof(basetype)]))
                  & bitpacker.max_overflow;
        } else
        {
          if constexpr (sizeof_unit == sizeof(basetype) + 2)
          {
            return ((integer << bitpacker.overflow_left_shift) |
                    (static_cast<uint64_t>(bytes[sizeof(basetype)]) << 8) |
                    (static_cast<uint64_t>(bytes[sizeof(basetype)+1])))
                    & bitpacker.max_overflow;
          } else
          {
            if constexpr (sizeof_unit == sizeof(basetype) + 3)
            {
              return ((integer << bitpacker.overflow_left_shift) |
                      (static_cast<uint64_t>(bytes[sizeof(basetype)]) << 16) |
                      (static_cast<uint64_t>(bytes[sizeof(basetype)+1]) << 8) |
                      (static_cast<uint64_t>(bytes[sizeof(basetype)+2])))
                      & bitpacker.max_overflow;
            } else
            {
              static_assert(sizeof_unit == sizeof(basetype) + 4);
              return ((integer << bitpacker.overflow_left_shift) |
                      (static_cast<uint64_t>(bytes[sizeof(basetype)]) << 24) |
                      (static_cast<uint64_t>(bytes[sizeof(basetype)+1]) << 16) |
                      (static_cast<uint64_t>(bytes[sizeof(basetype)+2]) << 8) |
                      (static_cast<uint64_t>(bytes[sizeof(basetype)+3])))
                      & bitpacker.max_overflow;
            }
          }
        }
#endif
      }
    }

#ifdef OWNCOPY_CONSTRUCTOR
    //It seems that this copy constructor is implicit and leaving it
    //leads to a compiler error or clang++. So we better exclude it here.
    BytesUnit& operator=(BytesUnit& other) noexcept
    {
      memcpy(&bytes[0],&other.bytes[0],sizeof_unit);
      return *this;
    }
#endif

    bool operator != (const BytesUnit& other) const noexcept
    {
      return memcmp(&bytes[0],&other.bytes[0],sizeof_unit) != 0;
    }

    bool operator == (const BytesUnit& other) const noexcept
    {
      return memcmp(&bytes[0],&other.bytes[0],sizeof_unit) == 0;
    }

    bool operator < (const BytesUnit& other) const noexcept
    {
      return memcmp(&bytes[0],&other.bytes[0],sizeof_unit) < 0;
    }

    [[nodiscard]] size_t sum(void) const noexcept
    {
      const basetype integer = *(reinterpret_cast<const basetype *>(bytes));

      if constexpr (sizeof_unit == sizeof(basetype))
      {
        return static_cast<size_t>(integer);
      } else
      {
        return static_cast<size_t>(integer) + bytes[sizeof(basetype)];
      }
    }

    void show(void) const noexcept
    {
      for (size_t j = 0; j < sizeof_unit; j++)
      {
        std::cout << std::right << std::setw(4)
                  << static_cast<int>(bytes[j])
                  << (j == sizeof_unit - 1 ? "\n" : " ");
      }
    }
};
#endif

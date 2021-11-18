#ifndef BYTES_UNIT_HPP
#define BYTES_UNIT_HPP
#include <cstdint>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <type_traits>
#include "utilities/bitpacker.hpp"

template<int sizeof_unit,int bit_groups>
class BytesUnit
{
  using basetype
    = typename std::conditional<sizeof_unit >= 8,
                                uint64_t,
                                typename std::conditional<sizeof_unit >= 4,
                                                          uint32_t,
                                                          uint16_t>::type>
                                ::type;
  private:
    uint8_t bytes[sizeof_unit];
  public:
    BytesUnit() {};
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
        if constexpr (sizeof_unit == sizeof(basetype) + 1)
        {
          bytes[sizeof(basetype)] = static_cast<uint8_t>(overflow_value);
        } else
        {
          bytes[sizeof(basetype)] = static_cast<uint8_t>(overflow_value >> 8);
          bytes[sizeof(basetype)+1] = static_cast<uint8_t>(overflow_value);
        }
      }
      memcpy(&bytes[0],reinterpret_cast<uint8_t *>(&integer),sizeof integer);
    }

    template<int idx>
    uint64_t decode_at(const GttlBitPacker<sizeof_unit,bit_groups>
                             &bitpacker) const noexcept
    {
      static_assert(idx >= 0 && idx < bit_groups);
      const basetype integer = *(reinterpret_cast<const basetype *>(bytes));

      if constexpr (sizeof_unit == sizeof(basetype) || idx < bit_groups - 1)
      {
        return static_cast<uint64_t>(integer >> bitpacker.shift_tab[idx]) &
               bitpacker.mask_tab[idx];
      } else /* sizeof_unit > sizeof(basetype) && idx == bit_groups - 1 */
      {
        if constexpr (sizeof_unit == sizeof(basetype) + 1)
        {
          return ((integer << bitpacker.overflow_left_shift) |
                  static_cast<uint64_t>(bytes[sizeof(basetype)]))
                  & bitpacker.max_overflow;
        } else
        {
          return ((integer << bitpacker.overflow_left_shift) |
                  (static_cast<uint64_t>(bytes[sizeof(basetype)]) << 8) |
                  (static_cast<uint64_t>(bytes[sizeof(basetype)+1])))
                  & bitpacker.max_overflow;
        }
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

    size_t sum(void) const noexcept
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

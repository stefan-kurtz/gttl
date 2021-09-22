#ifndef BYTES_UNIT_HPP
#define BYTES_UNIT_HPP
#include <cstdint>
#include <cstring>
#include <iostream>
#include <iomanip>
#include "utilities/bitpacker.hpp"

template<int sizeof_unit,int bit_groups>
class BytesUnit
{
  private:
    uint8_t bytes[sizeof_unit];
  public:
    BytesUnit() {};
    BytesUnit(const GttlBitPacker<sizeof_unit,bit_groups> &bitpacker,
              const std::array<uint64_t, bit_groups> &to_be_encoded)
    {
      static_assert(sizeof_unit >= sizeof(uint64_t));
      static_assert(sizeof *this == sizeof_unit);
      uint64_t integer = 0;
      if constexpr (sizeof_unit == 8)
      {
        for (int idx = 0; idx < bit_groups; idx++)
        {
          assert(to_be_encoded[idx] <= bitpacker.mask_tab[idx]);
          integer |= (to_be_encoded[idx] << bitpacker.shift_tab[idx]);
        }
      } else
      {
        static_assert (sizeof_unit > 8);
        for (int idx = 0; idx < bit_groups - 1; idx++)
        {
          assert(to_be_encoded[idx] <= bitpacker.mask_tab[idx]);
          integer |= (to_be_encoded[idx] << bitpacker.shift_tab[idx]);
        }
        const uint64_t overflow_value = to_be_encoded[bit_groups-1];
        assert(overflow_value <= bitpacker.max_overflow);
        integer |= (overflow_value >> bitpacker.overflow_left_shift);
        bytes[8] = static_cast<uint8_t>(overflow_value);
      }
      memcpy(&bytes[0],reinterpret_cast<uint8_t *>(&integer),8);
    }

    template<int idx>
    uint64_t decode_at(const GttlBitPacker<sizeof_unit,bit_groups> &bitpacker)
                       const noexcept
    {
      static_assert(idx >= 0 && idx < bit_groups);
      const uint64_t integer = *(reinterpret_cast<const uint64_t *>(bytes));

      if constexpr (sizeof_unit == 8 || idx < bit_groups - 1)
      {
        return (integer >> bitpacker.shift_tab[idx]) &
                bitpacker.mask_tab[idx];
      } else /* sizeof_unit == 9 && idx == bit_groups - 1 */
      {
        return ((integer << bitpacker.overflow_left_shift) |
                static_cast<uint64_t>(bytes[8])) & bitpacker.max_overflow;
      }
    }

    uint64_t decode_at0(const GttlBitPacker<sizeof_unit,bit_groups> &bitpacker)
                        const noexcept
    {
      return decode_at<0>(bitpacker);
    }

    uint64_t decode_at1(const GttlBitPacker<sizeof_unit,bit_groups> &bitpacker)
                        const noexcept
    {
      return decode_at<1>(bitpacker);
    }

    uint64_t decode_at2(const GttlBitPacker<sizeof_unit,bit_groups> &bitpacker)
                        const noexcept
    {
      return decode_at<2>(bitpacker);
    }

    uint64_t decode_at3(const GttlBitPacker<sizeof_unit,bit_groups> &bitpacker)
                        const noexcept
    {
      return decode_at<3>(bitpacker);
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
      const uint64_t integer = *(reinterpret_cast<const uint64_t *>(bytes));

      if constexpr (sizeof_unit == 8)
      {
        return static_cast<size_t>(integer);
      } else
      {
        return static_cast<size_t>(integer) + bytes[8];
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

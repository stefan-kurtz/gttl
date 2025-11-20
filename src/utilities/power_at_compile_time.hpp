#ifndef POWER_AT_COMPILE_TIME_HPP
#define POWER_AT_COMPILE_TIME_HPP
#include <cstdint>

template <typename T>
consteval T power_at_compile_time(T num, uint32_t exponent)
{
    return (exponent >= sizeof(uint32_t) * 8)
             ? 0 : exponent == 0 ? static_cast<T>(1)
                                 : num *
                                   power_at_compile_time(num, exponent - 1);
}
#endif

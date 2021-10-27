#ifndef CONSTEXPR_FOR_HPP
#define CONSTEXPR_FOR_HPP
#include <type_traits>
template <auto start, auto end, auto inc, class F>
constexpr void constexpr_for(F&& f)
{
  if constexpr (start < end)
  {
    f(std::integral_constant<decltype(start), start>());
    constexpr_for<start + inc, end, inc>(f);
  }
}
#endif

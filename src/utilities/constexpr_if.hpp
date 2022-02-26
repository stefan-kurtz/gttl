#ifndef CONSTEXPR_IF_HPP
#define CONSTEXPR_IF_HPP
#include <type_traits>
#include <utility>
#include <cstdbool>

/* from https://stackoverflow.com/questions/41011900/
                equivalent-ternary-operator-for-constexpr-if */

template <bool cond_v, typename Then, typename OrElse>
decltype(auto) constexpr_if(Then&& then, OrElse&& or_else)
{
  if constexpr (cond_v)
  {
    return std::forward<Then>(then);
  } else
  {
    return std::forward<OrElse>(or_else);
  }
}
#endif

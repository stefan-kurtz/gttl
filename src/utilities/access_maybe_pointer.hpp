#ifndef ACCESS_MAYBE_POINTER_HPP
#define ACCESS_MAYBE_POINTER_HPP
#include <type_traits>
#include <utility>

template <typename T>
decltype(auto) access_maybe_pointer(T &&t)
{
  if constexpr(std::is_pointer_v<std::remove_reference_t<T>>
               || requires { t.operator->(); })
  {
    return *t;
  } else {
    return std::forward<T>(t);
  }
}

#endif // ACCESS_MAYBE_POINTER_HPP

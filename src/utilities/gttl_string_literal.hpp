#ifndef GTTL_STRING_LITERAL_HPP
#define GTTL_STRING_LITERAL_HPP
#include <cstddef>
#include <algorithm>

template<size_t length_of_literal>
struct GttlStringLiteral
{
  /*
   Literal class type that wraps a constant expression string.
   Uses implicit conversion to allow templates to *seemingly* accept
   constant strings. Adapted from
   https://ctrpeach.io/posts/cpp20-string-literal-template-parameters/
  */
  char value[length_of_literal];
  constexpr GttlStringLiteral(const char (&str)[length_of_literal])
  {
    std::copy_n(str, length_of_literal, value);
  }
};
#endif

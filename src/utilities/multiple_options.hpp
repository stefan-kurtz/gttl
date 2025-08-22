#ifndef MULTIPLE_OPTIONS_HPP
#define MULTIPLE_OPTIONS_HPP

#include "utilities/find_lit_string.hpp"
#include "utilities/string_tokenizer.hpp"
#include "utilities/string_values_join.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

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


template<const GttlLitStringInitializerList &option_args>
class GttlMultipleOptions
{
  [[nodiscard]] uint64_t the_bit(int shift) const noexcept
  {
    return uint64_t{1} << shift;
  }

  uint64_t flags;
  const std::string option_name;

  public:
  explicit GttlMultipleOptions(const char* _option_name)
    : flags(0)
    , option_name(_option_name)
  { }

  [[nodiscard]] std::string help_string() const noexcept
  {
    return string_values_join(", ", option_args.begin(), option_args.end());
  }

  void set_flags(const std::string &argstring)
  {
    const StringTokenizer st(argstring);
    std::vector<std::string> option_args_string_vector(option_args.begin(),
                                                       option_args.end());

    for (auto && arg : st.token_list_get())
    {
      auto found = std::ranges::find(option_args_string_vector, arg);
      if (found == option_args_string_vector.end())
      {
        throw std::invalid_argument(std::string("illegal argument \"") +
                                    arg +
                                    std::string("\"") +
                                    option_name +
                                    std::string(", possible values: ") +
                                    help_string());
      }
      const size_t shift = std::distance(option_args_string_vector.begin(),
                                         found);
      this->flags |= the_bit(shift);
    }
  }

  template<GttlStringLiteral key>
  [[nodiscard]] bool is_set() const noexcept
  {
    constexpr const int shift
      = gttl_find_lit_string_at_compile_time(option_args, key.value);
    static_assert(shift < option_args.size());
    return flags & the_bit(shift);
  }

  [[nodiscard]] std::string to_string() const noexcept
  {
    std::vector<std::string> fields;
    int shift = 0;
    for (const auto *arg : option_args)
    {
      if (flags & the_bit(shift))
      {
        fields.emplace_back(arg);
      }
      ++shift;
    }
    return string_values_join(", ", fields.begin(), fields.end());
  }
};
#endif // MULTIPLE_OPTIONS_HPP

#ifndef BIT_SEQUENCE_HPP
#define BIT_SEQUENCE_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include <climits>

template<typename Basetype>
static inline std::string bit_sequence2string(Basetype bs)
{
  const unsigned int bits = CHAR_BIT * sizeof(Basetype);
  Basetype mask;
  std::string buffer{};

  buffer.reserve(bits);
  for (mask = (static_cast<Basetype>(1) << (bits - 1)); mask > 0; mask >>= 1)
  {
    buffer += (bs & mask) ? '1' : '0';
  }
  return buffer;
}

template<typename Basetype>
static inline std::string bit_sequence2string(Basetype bs,size_t group_size)
{
  static constexpr const unsigned int bits = CHAR_BIT * sizeof(Basetype);
  Basetype mask;
  std::string buffer{};
  buffer.reserve(bits + bits/group_size - 1);
  size_t gs = 0;
  for (mask = (static_cast<Basetype>(1) << (bits - 1)); mask > 0; mask >>= 1)
  {
    if (gs < group_size)
    {
      gs++;
    } else
    {
      buffer += ' ';
      gs = 1;
    }
    buffer += (bs & mask) ? '1' : '0';
  }
  return buffer;
}

template<typename Basetype>
static inline std::vector<std::string>
   bit_sequence2string_vector(const Basetype *encoding,size_t num_bytes)
{
  std::vector<std::string> s_vec{};
  s_vec.reserve(num_bytes);
  for (size_t idx = 0; idx < num_bytes; idx++)
  {
    auto bs = bit_sequence2string<uint8_t>(encoding[idx]);
    s_vec.push_back(bs);
  }
  return s_vec;
}
#endif

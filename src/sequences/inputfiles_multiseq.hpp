#ifndef INPUTFILES_MULTISEQ_HPP
#define INPUTFILES_MULTISEQ_HPP
#include <cstdint>
#include <exception>
#include <vector>
#include <string>
#include "sequences/gttl_multiseq.hpp"
#include "sequences/complement_plain.hpp"

static inline GttlMultiseq *
  gttl_inputfiles_multiseq(const std::vector<std::string> &inputfiles,
                           bool store_header,
                           bool store_sequence,
                           uint8_t padding_char,
                           bool reverse_complement_option)
{
  GttlMultiseq *multiseq = nullptr;
  try
  {
    if (reverse_complement_option)
    {
      multiseq = multiseq_with_reverse_complement<complement_plain>
                                                 (inputfiles,
                                                  store_header,
                                                  store_sequence,
                                                  padding_char);
    } else
    {
      multiseq = new GttlMultiseq(inputfiles,
                                  store_header,
                                  store_sequence,
                                  padding_char);
    }
  }
  catch (const std::exception &err)
  {
    delete multiseq;
    throw;
  }
  return multiseq;
}
#endif

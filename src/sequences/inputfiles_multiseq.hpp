#ifndef INPUTFILES_MULTISEQ_HPP
#define INPUTFILES_MULTISEQ_HPP
#include <cstdint>
#include <exception>
#include <vector>
#include <string>
#include "sequences/gttl_multiseq.hpp"

static inline GttlMultiseq *
  gttl_inputfiles_multiseq(const std::vector<std::string> &inputfiles,
                           bool store_header,
                           bool store_sequence,
                           uint8_t padding_char,
                           bool reverse_complement_option)
{
  // NOLINTNEXTLINE(misc-const-correctness)
  GttlMultiseq *multiseq = nullptr;
  try
  {
     multiseq = new GttlMultiseq(inputfiles, /* CONSTRUCTOR */
                                 store_header,
                                 store_sequence,
                                 padding_char,
                                 reverse_complement_option);
  }
  catch (const std::exception &err)
  {
    delete multiseq;
    throw;
  }
  return multiseq;
}
#endif

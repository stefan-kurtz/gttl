#ifndef MULTISEQ_PAIR_HPP
#define MULTISEQ_PAIR_HPP
#include <cstdint>
#include <cstring>
#include <cassert>
#include <exception>
#include <stdexcept>
#include <string>
#include <tuple>
#include <format>
#include "sequences/gttl_multiseq.hpp"
#include "sequences/guess_if_protein_seq.hpp"

static inline std::tuple<GttlMultiseq *,GttlMultiseq *,bool>
  create_multiseq_pair(const char *dbfile,const char*queryfile)
{
  GttlMultiseq *db_multiseq = nullptr;
  // NOLINTNEXTLINE(misc-const-correctness)
  GttlMultiseq *query_multiseq = nullptr;
  bool is_dna_alphabet = false;
  static constexpr const bool store_header = true;
  static constexpr const bool store_sequence = true;
  constexpr const bool with_reverse_complement = false;
  assert(dbfile != nullptr);
  try
  {
    constexpr const uint8_t padding_char = UINT8_MAX;
    db_multiseq = new GttlMultiseq(std::string(dbfile), /* CONSTRUCTOR*/
                                   store_header,
                                   store_sequence,
                                   padding_char,
                                   with_reverse_complement);
  }
  catch (const std::exception &err)
  {
    throw;
  }
  assert(queryfile != nullptr);
  if (dbfile == queryfile || strcmp(dbfile,queryfile) == 0)
  {
    query_multiseq = db_multiseq;
  } else
  {
    try
    {
      constexpr const uint8_t padding_char = UINT8_MAX - 1;
      query_multiseq = new GttlMultiseq(std::string(queryfile), /* CONSTRUCTOR*/
                                        store_header,
                                        store_sequence,
                                        padding_char,
                                        with_reverse_complement);
    }
    catch (const std::exception &err)
    {
      if (db_multiseq != query_multiseq)
      {
        delete query_multiseq;
      }
      delete db_multiseq;
      throw;
    }
  }
  if (guess_if_protein_multiseq(db_multiseq))
  {
    if (query_multiseq != db_multiseq &&
        !guess_if_protein_multiseq(query_multiseq))
    {
      if(db_multiseq != query_multiseq)
      {
        delete query_multiseq;
      }
      delete db_multiseq;
      throw std::runtime_error(
              std::format(": incompatible files: file \"{}\" contains protein "
                          "sequences, but file \"{}\" does not",
                          dbfile,
                          queryfile));
    }
  } else
  {
    is_dna_alphabet = true;
    if (query_multiseq != db_multiseq &&
        guess_if_protein_multiseq(query_multiseq))
    {
      if(db_multiseq != query_multiseq)
      {
        delete query_multiseq;
      }
      delete db_multiseq;
      throw std::runtime_error(
              std::format(": incompatible files: file \"{}\" does not contain "
                          "protein sequences, but file \"{}\" does",
                          dbfile,
                          queryfile));
    }
  }
  return {db_multiseq, query_multiseq, is_dna_alphabet};
}
#endif

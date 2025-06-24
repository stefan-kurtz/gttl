#ifndef MULTISEQ_PAIR_HPP
#define MULTISEQ_PAIR_HPP
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <exception>
#include <stdexcept>
#include <tuple>
#include "sequences/gttl_multiseq.hpp"
#include "sequences/guess_if_protein_seq.hpp"
#include "utilities/str_format.hpp"

static inline std::tuple<GttlMultiseq *,GttlMultiseq *,bool>
  create_multiseq_pair(const char *dbfile,const char*queryfile)
{
  GttlMultiseq *db_multiseq = NULL;
  GttlMultiseq *query_multiseq = NULL;
  bool is_dna_alphabet = false;
  static constexpr const bool store_header = true;
  static constexpr const bool store_sequence = true;
  assert(dbfile != NULL);
  try
  {
    db_multiseq = new GttlMultiseq(dbfile,store_header,store_sequence,
                                   UINT8_MAX);
  }
  catch (const std::exception &err)
  {
    throw;
  }
  assert(queryfile != NULL);
  if (dbfile == queryfile || strcmp(dbfile,queryfile) == 0)
  {
    query_multiseq = db_multiseq;
  } else
  {
    try
    {
      query_multiseq = new GttlMultiseq(queryfile,store_header,store_sequence,
                                        UINT8_MAX-1);
    }
    catch (const std::exception &err)
    {
      if(db_multiseq != query_multiseq) delete query_multiseq;
      delete db_multiseq;
      throw;
    }
  }
  if (guess_if_protein_multiseq(db_multiseq))
  {
    if (query_multiseq != db_multiseq &&
        !guess_if_protein_multiseq(query_multiseq))
    {
      if(db_multiseq != query_multiseq) delete query_multiseq;
      delete db_multiseq;
      const StrFormat msg(": incompatible files: file \"%s\" contains protein "
                          "sequences, but file \"%s\" does not",
                          dbfile,
                          queryfile);
      throw std::runtime_error(msg.str());
    }
  } else
  {
    is_dna_alphabet = true;
    if (query_multiseq != db_multiseq &&
        guess_if_protein_multiseq(query_multiseq))
    {
      if(db_multiseq != query_multiseq) delete query_multiseq;
      delete db_multiseq;
      const StrFormat msg(": incompatible files: file \"%s\" does not contain "
                          "protein sequences, but file \"%s\" does",
                          dbfile,
                          queryfile);
      throw std::runtime_error(msg.str());
    }
  }
  return std::tuple<GttlMultiseq *,GttlMultiseq *,bool>
                   (db_multiseq,query_multiseq,is_dna_alphabet);
}
#endif

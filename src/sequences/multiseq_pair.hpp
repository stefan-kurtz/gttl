#ifndef MULTISEQ_PAIR_HPP
#define MULTISEQ_PAIR_HPP
#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstring>
#include <cassert>
#include <tuple>
#include "sequences/gttl_multiseq.hpp"
#include "sequences/guess_if_protein_seq.hpp"

static inline std::tuple<GttlMultiseq *,GttlMultiseq *,bool>
  create_multiseq_pair(const char *dbfile,const char*queryfile)
{
  GttlMultiseq *db_multiseq = NULL, *query_multiseq = NULL;
  bool dna_alphabet = false;
  static constexpr const bool store_sequences = true;
  assert(dbfile != NULL);
  try
  {
    db_multiseq = new GttlMultiseq(dbfile,store_sequences,UINT8_MAX);
  }
  catch (std::string &msg)
  {
    throw msg;
  }
  assert(queryfile != NULL);
  if (dbfile == queryfile || strcmp(dbfile,queryfile) == 0)
  {
    query_multiseq = db_multiseq;
  } else
  {
    try
    {
      query_multiseq = new GttlMultiseq(queryfile,store_sequences,UINT8_MAX-1);
    }
    catch (std::string &msg)
    {
      throw msg;
    }
  }
  if (guess_if_protein_multiseq(db_multiseq))
  {
    if (query_multiseq != db_multiseq &&
        !guess_if_protein_multiseq(query_multiseq))
    {
      StrFormat msg(": incompatible files: file \"%s\" contains protein "
                    "sequences, but file \"%s\" does not",
                    dbfile,queryfile);
      throw msg.str();
    }
  } else
  {
    dna_alphabet = true;
    if (query_multiseq != db_multiseq &&
        guess_if_protein_multiseq(query_multiseq))
    {
      StrFormat msg(": incompatible files: file \"%s\" does not contain "
                    "protein sequences, but file \"%s\" does",
                    dbfile,queryfile);
      throw msg.str();
    }
  }
  return std::tuple<GttlMultiseq *,GttlMultiseq *,bool>
                   (db_multiseq,query_multiseq,dna_alphabet);
}

static inline std::tuple<GttlMultiseq *,GttlMultiseq *,bool>
  create_multiseq_pair(const std::string &dbfile,const std::string &queryfile)
{
  return create_multiseq_pair(dbfile.c_str(),queryfile.c_str());
}
#endif

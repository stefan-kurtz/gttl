#ifndef SUCCINCT_PLCP_TABLE_HPP
#define SUCCINCT_PLCP_TABLE_HPP

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include "utilities/gttl_binary_read.hpp"
#include "indexes/succinct_bitvector.hpp"

template <typename SuftabBaseType>
class SuccinctPlcpTable
{
  class Iterator
  {
    const BinaryFileReader<SuftabBaseType> suftab_reader;
    BinaryFileReader<SuftabBaseType>::Iterator suftab_iter;
    BinaryFileReader<SuftabBaseType>::Iterator suftab_end;
    size_t totallength, idx;
    uint32_t current_lcpvalue;
    const SuccinctBitvector *succinct_plcp_table;
    public:
    Iterator(size_t _totallength,
             size_t _idx,
             const std::string &suftab_inputfile,
             const SuccinctBitvector *_succinct_plcp_table)
      : suftab_reader(suftab_inputfile)
      , suftab_iter(suftab_reader.begin())
      , totallength(_totallength)
      , idx(_idx)
      , current_lcpvalue(0)
      , succinct_plcp_table(_succinct_plcp_table)
    {
      next();
    }
    uint32_t operator*(void) const
    {
      return current_lcpvalue;
    }
    Iterator& operator++(void)
    {
      next();
      return *this;
    }
    void next(void) {
      ++suftab_iter;
      ++idx;
      if (idx < totallength)
      {
        const SuftabBaseType current_suftab = *suftab_iter;
        const uint32_t suffix               = current_suftab + 1;
        const size_t select_1 = succinct_plcp_table->get_select(suffix, 1);
        const size_t lcp      = succinct_plcp_table->get_rank(select_1, 0) + 1
                                     - suffix;
        current_lcpvalue = lcp;
      } else
      {
        //printf("%zu\n", idx);
        current_lcpvalue = 0;
      }
    }
    bool operator != (const Iterator& other) const
    {
      return idx != other.idx;
    }
  };
  const std::string suftab_inputfile;
  const std::string lcp_inputfile;
  SuccinctBitvector succinct_plcp_table;
  size_t totallength;
  public:
  SuccinctPlcpTable(const std::string &indexname)
   : suftab_inputfile(indexname + ".suf")
   , lcp_inputfile(indexname + ".lls")
   , succinct_plcp_table(lcp_inputfile.c_str())
   , totallength(succinct_plcp_table.get_rank(succinct_plcp_table.get_length(),
                                              1))
  { }
  ~SuccinctPlcpTable(void)
  { }
  [[nodiscard]] Iterator begin(void) const
  {
    return Iterator(totallength, 0, suftab_inputfile, &succinct_plcp_table);
  }
  [[nodiscard]] Iterator end(void) const
  {
    return Iterator(0, totallength, suftab_inputfile, nullptr);
  }
};
#endif

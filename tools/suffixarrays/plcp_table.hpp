#ifndef PLCP_TABLE_HPP
#define PLCP_TABLE_HPP

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <string>
#include "utilities/memory_tracker.hpp"
#include "utilities/gttl_binary_read.hpp"

template <typename SuftabBaseType>
class PlcpTable
{
  class Iterator
  {
    const BinaryFileReader<SuftabBaseType> suftab_reader;
    BinaryFileReader<SuftabBaseType>::Iterator suftab_iter;
    size_t totallength, idx;
    uint32_t current_lcpvalue;
    const uint32_t *plcp_table;
    public:
    Iterator(size_t _totallength,
             size_t _idx,
             const std::string &suftab_inputfile,
             const uint32_t *_plcp_table)
      : suftab_reader(suftab_inputfile)
      , suftab_iter(suftab_reader.begin())
      , totallength(_totallength)
      , idx(_idx)
      , current_lcpvalue(0)
      , plcp_table(_plcp_table)
    { }
    uint32_t operator*(void) const
    {
      return current_lcpvalue;
    }
    Iterator& operator++(void)
    {
      ++suftab_iter;
      ++idx;
      if (idx < totallength)
      {
        const SuftabBaseType current_suftab = *suftab_iter;
        current_lcpvalue = plcp_table[current_suftab];
      } else
      {
        current_lcpvalue = 0;
      }
      return *this;
    }
    bool operator != (const Iterator& other) const
    {
      return idx != other.idx;
    }
  };
  GttlMemoryTracker *memory_tracker;
  const std::string suftab_inputfile;
  size_t totallength;
  SuftabBaseType *phi_table;
  uint32_t *plcp_table;
  public:
  PlcpTable(GttlMemoryTracker *_memory_tracker,
            const uint8_t *sequence,
            size_t _totallength,
            const std::string &indexname,
            size_t alphasize)
   : memory_tracker(_memory_tracker)
   , suftab_inputfile(indexname + ".suf")
   , totallength(_totallength)
   , phi_table(new SuftabBaseType [totallength + 1])
     /* reuse space for plcp_table */
   , plcp_table(reinterpret_cast<uint32_t *>(phi_table))
  {
    memory_tracker->track(phi_table, __FILE__, __LINE__,
                          (totallength + 1) * sizeof *phi_table);
    const BinaryFileReader<SuftabBaseType> snd_suftab_reader(suftab_inputfile);
    auto suftab_iter = snd_suftab_reader.begin();
    SuftabBaseType previous_suftab = *suftab_iter;
    phi_table[previous_suftab] = static_cast<SuftabBaseType>(totallength);
    for(++suftab_iter; suftab_iter != snd_suftab_reader.end(); ++suftab_iter)
    {
      const SuftabBaseType current_suftab = *suftab_iter;
      assert(current_suftab <= totallength);
      /* random access to phi_table */
      phi_table[current_suftab] = previous_suftab;
      previous_suftab = current_suftab;
    }
    uint32_t lcpvalue = 0;
    for (size_t pos = 0; pos < totallength; pos++)
    {
      /* sequential access to phi_table which is in memory */
      const size_t previousstart = phi_table[pos];
      /* sequential access to plcp_table */
      while (pos + lcpvalue < totallength &&
             previousstart + lcpvalue < totallength)
      {
        const uint8_t cc0 = sequence[pos + lcpvalue];
        const uint8_t cc1 = sequence[previousstart + lcpvalue];
        if (cc0 != cc1 or cc0 >= static_cast<uint8_t>(alphasize))
        {
          break;
        }
        lcpvalue++;
      }
      plcp_table[pos] = lcpvalue;
      lcpvalue -= (lcpvalue > 0);
    }
  }
  uint32_t plcp_value_get(size_t index) const
  {
    return plcp_table[index];
  }
  size_t get_total_length() const
  {
    return totallength;
  }
  ~PlcpTable(void)
  {
    memory_tracker->untrack(phi_table, __FILE__, __LINE__);
    delete[] phi_table;
  }
  Iterator begin(void) const
  {
    return Iterator(totallength, 0, suftab_inputfile, plcp_table);
  }
  Iterator end(void) const
  {
    return Iterator(0, totallength + 1, suftab_inputfile, nullptr);
  }
};
#endif

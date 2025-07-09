/*
  Copyright (c) 2022-2025 Stefan Kurtz <stefan.kurtz@uni-hamburg.de>
  Copyright (c) 2022-2025 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/
#ifndef GTTL_SUFFIXARRAY_HPP
#define GTTL_SUFFIXARRAY_HPP

#include <cassert>
#include <cstddef>
#include <string>
#include <fstream>
#include <ios>
#include <cstdio>
#include <vector>
#include <cstdint>
#include <cinttypes>
#include <tuple>
#include <algorithm>

#include "utilities/read_vector.hpp"
#include "utilities/gttl_mmap.hpp"

class LCPtable
{
  private:
  struct Iterator
  {
    private:
    const std::vector<uint8_t> &small_lcptab;
    const std::vector<uint16_t> &ll2tab;
    const std::vector<uint32_t> &ll4tab;
    size_t current_idx, end_idx, ll2_idx, ll4_idx;
    uint32_t current_lcpvalue;
    uint32_t current_lcpvalue_get(size_t idx)
    {
      assert(idx < small_lcptab.size());
      if (small_lcptab[idx] != UINT8_MAX)
      {
        return static_cast<uint32_t>(small_lcptab[idx]);
      }
      assert(ll2_idx < ll2tab.size());
      if (ll2tab[ll2_idx] < UINT16_MAX)
      {
        return static_cast<uint32_t>(ll2tab[ll2_idx++]);
      }
      ll2_idx++;
      assert(ll4_idx < ll4tab.size());
      return ll4tab[ll4_idx++];
    }
    public:
    Iterator(const std::vector<uint8_t> &_small_lcptab,
             const std::vector<uint16_t> &_ll2tab,
             const std::vector<uint32_t> &_ll4tab)
      : small_lcptab(_small_lcptab)
      , ll2tab(_ll2tab)
      , ll4tab(_ll4tab)
      , current_idx(1)
      , end_idx(_small_lcptab.size())
      , ll2_idx(0)
      , ll4_idx(0)
      , current_lcpvalue(current_lcpvalue_get(1))
    { }
    Iterator& operator++()
    {
      current_idx++;
      assert(current_idx < end_idx);
      current_lcpvalue = current_lcpvalue_get(current_idx);
      return *this;
    }
    uint32_t operator*(void) const noexcept
    {
      return current_lcpvalue;
    }
    bool operator != (const Iterator& other) const noexcept
    {
      return current_idx != other.end_idx;
    }
  };
  std::vector<uint8_t> small_lcptab;
  std::vector<uint16_t> ll2tab;
  std::vector<uint32_t> ll4tab;
public:
  LCPtable(const std::string &infile_base)
    : small_lcptab(gttl_read_vector<uint8_t>(infile_base + ".lcp"))
    , ll2tab(gttl_read_vector<uint16_t>(infile_base + ".ll2"))
    , ll4tab(gttl_read_vector<uint32_t>(infile_base + ".ll4"))
  { }

  [[nodiscard]] Iterator begin(void) const noexcept
  {
    return Iterator(small_lcptab, ll2tab, ll4tab);
  }
  [[nodiscard]] Iterator end(void) const noexcept
  {
    return Iterator(small_lcptab, ll2tab, ll4tab);
  }
};

enum Suffixarrayfiles : uint8_t
{
  LCPTAB_file,
  SUFTAB_file,
  BU_SUFTAB_file,
  MMAP_BU_SUFTAB_file,
  TIS_file
};

class GttlSuffixArray
{
  using SuftabBaseType = uint32_t;
  /* Adjust the following value, when adding additional integer keys */
  static constexpr const int num_integer_keys = 6;
  const std::vector<std::string> keys = {"reverse_complement",
                                         "nonspecial_suffixes",
                                         "sequences_number",
                                         "sequences_number_bits",
                                         "sequences_length_bits",
                                         "sizeof_suftab_entry",
                                         "inputfile"};
  private:
  std::vector<SuftabBaseType> suftab_abspos;
  std::vector<uint8_t> suftab_bytes;
  Gttlmmap<uint8_t> *gttl_mmap_suftab_bytes;

  std::vector<uint8_t> tistab;
  LCPtable *lcptable;
  size_t int_values[num_integer_keys];
  bool int_values_set[num_integer_keys] = {false};
  std::vector<std::string> inputfiles;

  [[nodiscard]] int key2index(const std::string &key) const noexcept
  {
    auto found = std::find(keys.begin(), keys.end(), key);
    if (found == keys.end())
    {
      return -1;
    }
    return static_cast<int>(found - keys.begin());
  }

  void read_in_prj_file(const std::string &prj_filename)
  {
    std::string line;
    size_t sep_pos = 0;
    std::ifstream in_file;

    in_file.open(prj_filename, std::ifstream::in);
    if (in_file.fail() )
    {
      throw std::ios_base::failure(std::string("file ").append(prj_filename)
                                   .append(": cannot open; possibly index "
                                           "needs to be created"));
    }
    /* number of keys with integer values */
    while (std::getline (in_file,line))
    {
      if((sep_pos = line.find('\t')) == std::string::npos)
      {
        throw std::ios_base::failure(std::string("file ").append(prj_filename)
                                     .append(": missing tabulator in line ")
                                     .append(line));
      }
      const std::string this_key = line.substr(0, sep_pos);
      const int idx = key2index(this_key);
      if (idx == -1)
      {
        throw std::ios_base::failure(std::string("file ").append(prj_filename)
                                     .append(": illegal key ").append(this_key)
                                     .append(" in line ").append(line));
      }
      if (idx < num_integer_keys)
      {
        int64_t read_int;
        if (sscanf(line.substr(sep_pos + 1).c_str(),"%" PRId64, &read_int) != 1
            || read_int < 0)
        {
          throw std::ios_base::failure(std::string("file ").append(prj_filename)
                                    .append(": value for key ")
                                    .append(this_key)
                                    .append(" must be a positive integer"));
        }
        int_values[idx] = static_cast<size_t>(read_int);
        int_values_set[idx] = true;
      } else
      {
        inputfiles.push_back(line.substr(sep_pos + 1));
      }
    }
    in_file.close();
    for (int idx = 0; idx < num_integer_keys; idx++)
    {
      if (!int_values_set[idx])
      {
        throw std::ios_base::failure(std::string("file ") + prj_filename +
                                     std::string(": missing line for integer "
                                                 "key \"") +
                                     keys[idx] + std::string("\""));
      }
    }
    if (inputfiles.empty())
    {
      throw std::ios_base::failure(std::string("file ") + prj_filename +
                                   std::string(": missing lines with key "
                                               "\"inputfile\""));
    }
  }
  public:
  GttlSuffixArray(const std::string &infile_base,
                  const std::vector<Suffixarrayfiles> &saf_vec)
    : gttl_mmap_suftab_bytes(nullptr)
    , lcptable(nullptr)
  {
    read_in_prj_file(infile_base + ".prj");
    for(const auto& value: saf_vec)
    {
      if (value == LCPTAB_file)
      {
        lcptable = new LCPtable(infile_base);
        continue;
      }
      if (value == SUFTAB_file)
      {
        suftab_abspos = gttl_read_vector<SuftabBaseType>(infile_base + ".suf");
        continue;
      }
      if (value == BU_SUFTAB_file)
      {
        suftab_bytes = gttl_read_vector<uint8_t>(infile_base + ".bsf");
        continue;
      }
      if (value == MMAP_BU_SUFTAB_file)
      {
        const std::string bu_inputfile{infile_base + ".bsf"};
        gttl_mmap_suftab_bytes = new Gttlmmap<uint8_t>(bu_inputfile.c_str());
        continue;
      }
      if (value == TIS_file)
      {
        tistab = gttl_read_vector<uint8_t>(infile_base + ".tis");
        continue;
      }
    }
  }
  ~GttlSuffixArray()
  {
    if (lcptable != nullptr)
    {
      delete lcptable;
    }
    if (gttl_mmap_suftab_bytes != nullptr)
    {
      delete gttl_mmap_suftab_bytes;
    }
  }
  [[nodiscard]] const std::vector<SuftabBaseType> &
  get_suftab_abspos() const noexcept
  {
    assert(not suftab_abspos.empty());
    return suftab_abspos;
  }
  [[nodiscard]] const std::vector<uint8_t> &get_suftab_bytes() const noexcept
  {
    assert(not suftab_bytes.empty());
    return suftab_bytes;
  }
  [[nodiscard]] const uint8_t *get_mmap_suftab_bytes(void) const noexcept
  {
    assert(gttl_mmap_suftab_bytes);
    return gttl_mmap_suftab_bytes->ptr();
  }
  [[nodiscard]] const std::vector<uint8_t> &get_tistab(void) const noexcept
  {
    assert(not tistab.empty());
    return tistab;
  }
  [[nodiscard]] const LCPtable &lcptable_get(void) const noexcept
  {
    assert(lcptable != nullptr);
    return *lcptable;
  }
  [[nodiscard]] bool with_reverse_complement(void) const noexcept
  {
    return int_values[key2index("reverse_complement")] == 1 ? true : false;
  }
  [[nodiscard]] size_t nonspecial_suffixes_get(void) const noexcept
  {
    return int_values[key2index("nonspecial_suffixes")];
  }
  [[nodiscard]] size_t sequences_number_get(void) const noexcept
  {
    return int_values[key2index("sequences_number")];
  }
  [[nodiscard]] int sequences_number_bits_get(void) const noexcept
  {
    return static_cast<int>(int_values[key2index("sequences_number_bits")]);
  }
  [[nodiscard]] int sequences_length_bits_get(void) const noexcept
  {
    return static_cast<int>(int_values[key2index("sequences_length_bits")]);
  }
  [[nodiscard]] int sizeof_suftab_entry(void) const noexcept
  {
    return static_cast<int>(int_values[key2index("sizeof_suftab_entry")]);
  }
  [[nodiscard]] const std::vector<std::string> &
  inputfiles_get(void) const noexcept
  {
    return inputfiles;
  }
  private:
  [[nodiscard]] SuftabBaseType lcp_interval_find_rightbound(uint8_t cc,
                                                          size_t offset,
                                                          SuftabBaseType left,
                                                          SuftabBaseType right)
                                                           const
  {
    while (left + 1 < right)
    {
      const SuftabBaseType mid = (left + right)/2;
      const size_t pos = offset + suftab_abspos[mid];
      const uint8_t midcc = tistab[pos];
      if (cc < midcc)
      {
        right = mid;
      } else
      {
        left = mid;
      }
    }
    return left;
  }
  [[nodiscard]] uint8_t sequence_access(size_t idx) const
  {
    return idx == tistab.size() ? UINT8_MAX : tistab[idx];
  }
  [[nodiscard]] std::tuple<bool,SuftabBaseType,SuftabBaseType>
    lcp_interval_find_child_intv(uint8_t cc,
                                 size_t offset,
                                 SuftabBaseType left,
                                 SuftabBaseType right) const
  {
    SuftabBaseType leftbound = left;
    size_t pos = offset + suftab_abspos[right];
    uint8_t leftcc;
    const uint8_t rightcc = sequence_access(pos);
    while (true)
    {
      pos = offset + suftab_abspos[leftbound];
      leftcc = sequence_access(pos);
      if (leftcc == rightcc)
      {
        break;
      }
      const SuftabBaseType rightbound
        = lcp_interval_find_rightbound(leftcc,offset,leftbound,right);
      if (leftcc == cc)
      {
        return std::make_tuple(true,leftbound,rightbound);
      }
      if (leftcc > cc)
      {
        return std::make_tuple(false,0,0);
      }
      leftbound = rightbound + 1;
    }
    if (leftcc == cc)
    {
      return std::make_tuple(true,leftbound,right);
    }
    return std::make_tuple(false,0,0);
  }
  public:
  auto find_maximal_prefix(const uint8_t *query, size_t querylen) const
  {
    const size_t nonspecial_suffixes = nonspecial_suffixes_get();

    assert(nonspecial_suffixes > 0 and nonspecial_suffixes <= UINT32_MAX);
    SuftabBaseType left = 0;
    SuftabBaseType right = nonspecial_suffixes - 1;
    size_t idx;
    for (idx = 0; idx < querylen; idx++)
    {
      auto result = lcp_interval_find_child_intv(query[idx],
                                                 idx,
                                                 left,
                                                 right);
      if (std::get<0>(result))
      {
        left = std::get<1>(result);
        right = std::get<2>(result);
      } else
      {
        break;
      }
    }
    return std::make_tuple(idx,suftab_abspos.begin() + left,
                               suftab_abspos.begin() + right + 1);
  }
  size_t find_maximal_prefix_length(const uint8_t *query, size_t querylen) const
  {
    return std::get<0>(find_maximal_prefix(query,querylen));
  }
};
#endif  // SUFFIXARRAY_HPP

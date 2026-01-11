/*
  Copyright (c) 2022 Stefan Kurtz <stefan.kurtz@uni-hamburg.de>
  Copyright (c) 2022 Center for Bioinformatics, University of Hamburg

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

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <format>
#include "utilities/runtime_class.hpp"
#include "indexes/gttl_suffixarray.hpp"
#include "indexes/succinct_bitvector.hpp"
#include "succinct_plcp_table.hpp"

int main(int argc,char *argv[])
{
  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " <indexname>\n";
    return EXIT_FAILURE;
  }
  const char *const indexname  = argv[1];
  const GttlSuffixArray *suffixarray = nullptr;
  RunTimeClass rt_overall{};
  bool haserr = false;
  try
  {
    suffixarray = new GttlSuffixArray(indexname,{LCPTAB_file, SUFTAB_file});
  }
  catch (const std::exception &err)
  {
    std::cerr << argv[0] << ": indexname \"" << indexname
              << "\": " << err.what() << '\n';
    haserr = true;
  }
  if (!haserr)
  {
    std::cout << "# read the succinct representaton from file " << indexname
              << ".lls\n";
    const std::string lls_filename{std::string(indexname) + ".lls"};
    const SuccinctBitvector succinctlcp(lls_filename);
    const SuccinctPlcpTable<uint32_t> succinctplcptable(indexname);
    auto succinctplcpiter = succinctplcptable.begin();
    const size_t nonspecial_suffixes = suffixarray->nonspecial_suffixes_get();
    const LCPtable &lcptable = suffixarray->lcptable_get();
    size_t idx = 0;
    size_t lcpvalue_sum = 0;
    for (const uint32_t lcpvalue : lcptable)
    {
       /* at idx */
      /* HAL: verify that the lcp value at index idx
         from the succinct representation equals lcpvalue */
      const uint32_t suffix = suffixarray->get_suftab_abspos().at(idx + 1) + 1;
      try
      {
        size_t lcp;
        if (suffix > (succinctlcp.length_get() + 1)  / 2) {
          lcp = 0;
        } else {
          const size_t select_1 = succinctlcp.get_select(suffix, true);
          lcp =  succinctlcp.get_rank(select_1, false) + 1 - suffix;
        }
        const uint32_t succinct_lcp = *succinctplcpiter;
        if (lcp != lcpvalue || succinct_lcp != lcpvalue)
        {
          fprintf(stderr,"lcpmismatch: %zu, %" PRIu32 ", %zu, %" PRIu32 "\n",
                         idx, lcpvalue, lcp, succinct_lcp);
          exit(EXIT_FAILURE);
        }
      }
      catch (const std::exception &err)
      {
        std::cerr << argv[0] << ": indexname \"" << indexname
                  << "\": " << err.what() << '\n';
        haserr = true;
        break;
      }

      ++succinctplcpiter;
      lcpvalue_sum += static_cast<size_t>(lcpvalue);
      idx++;
      if (idx == nonspecial_suffixes)
      {
        /* all following lcp-values are 0 as they are all suffixes
           beginning with special symbol which is different from all other
           symbols */
        break;
      }
    }
    std::cout << "# reverse_complement\t"
              << (suffixarray->with_reverse_complement() ? "true" : "false")
              << '\n';
    std::cout << "# average lcp value\t"
              << (static_cast<double>(lcpvalue_sum) / nonspecial_suffixes)
              << '\n';
  }
  delete suffixarray;
  if (!haserr)
  {
    rt_overall.show(std::format("lcp_check\t{}", indexname));
  }
  return haserr ? EXIT_FAILURE : EXIT_SUCCESS;
}

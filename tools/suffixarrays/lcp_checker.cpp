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
#include "utilities/runtime_class.hpp"
#include "indexes/gttl_suffixarray.hpp"
#include "indexes/succinct_bitvector.hpp"
#include "succinct_plcp_table.hpp"
#include "utilities/str_format.hpp"

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
    /* HAL: add your code here */
    std::cout << "read the succinct representaton from file " << indexname
              << ".lls\n";
  }
  const std::string filename          = std::string(indexname) + ".lls";
  const SuccinctBitvector succinctlcp = SuccinctBitvector(filename.c_str());

  const std::string indexstring(indexname);
  const SuccinctPlcpTable<uint32_t> succinctplcptable(indexstring);
  auto succinctplcpiter = succinctplcptable.begin();
  if (!haserr)
  {
    const size_t nonspecial_suffixes = suffixarray->nonspecial_suffixes_get();
    const LCPtable &lcptable = suffixarray->lcptable_get();
    size_t idx = 0;
    size_t lcpvalue_sum = 0;
    for (auto it = lcptable.begin(); it != lcptable.end(); ++it)
    {
      const uint32_t lcpvalue = *it; /* at idx */
      /* HAL: verify that the lcp value at index idx
         from the succinct representation equals lcpvalue */

      const uint32_t suffix = suffixarray->get_suftab_abspos().at(idx + 1) + 1;
      const size_t select_1 = succinctlcp.get_select(suffix, true);
      const size_t lcp      = succinctlcp.get_rank(select_1, false)
                              + 1 - suffix;


      const uint32_t succinct_lcp = *succinctplcpiter;
      if (lcp != lcpvalue || succinct_lcp != lcpvalue) {
        printf("lcpmismatch: %zu, %d, %zu, %d\n",
               idx, lcpvalue, lcp, succinct_lcp);
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

    // printf("%zu\n", idx);
    // for (; succinctplcpiter != succinctplcptable.end(); ++succinctplcpiter) {
    //   printf("additional %d\n", *succinctplcpiter);
    // }
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
    const StrFormat msg("sa_reader\t%s", indexname);
    rt_overall.show(msg.str());
  }
  return haserr ? EXIT_FAILURE : EXIT_SUCCESS;
}

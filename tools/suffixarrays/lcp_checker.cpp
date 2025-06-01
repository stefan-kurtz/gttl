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
#include <cstdlib>
#include <iostream>
#include "utilities/runtime_class.hpp"
#include "indexes/gttl_suffixarray.hpp"

int main(int argc,char *argv[])
{
  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " <indexname>"
              << std::endl;
    return EXIT_FAILURE;
  }
  const char *indexname = argv[1];
  GttlSuffixArray *suffixarray = nullptr;
  RunTimeClass rt_overall{};
  bool haserr = false;
  try
  {
    suffixarray = new GttlSuffixArray(indexname,{LCPTAB_file});
  }
  catch (const std::runtime_error &err)
  {
    std::cerr << argv[0] << ": indexname \"" << indexname << "\": "
              << err.what() << std::endl;
    haserr = true;
  }
  if (!haserr)
  {
    /* HAL: add your code here */
    std::cout << "read the succinct representaton from file " << indexname
              << ".lls" << std::endl;
  }
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
              << std::endl;
    std::cout << "# average lcp value\t" << (static_cast<double>(lcpvalue_sum)/
                                             nonspecial_suffixes)
              << std::endl;
  }
  delete suffixarray;
  if (!haserr)
  {
    StrFormat msg("sa_reader\t%s",indexname);
    rt_overall.show(msg.str());
  }
  return haserr ? EXIT_FAILURE : EXIT_SUCCESS;
}

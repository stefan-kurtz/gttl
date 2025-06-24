/*
  Copyright (c) 2021-2022 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2021-2022 Center for Bioinformatics, University of Hamburg

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
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include "sequences/guess_if_protein_seq.hpp"

int main(int argc,char *argv[])
{
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " <inputfile1> [inputfile2...]"
              << '\n';
    return EXIT_FAILURE;
  }
  std::vector<std::string> inputfiles{};
  for (int idx = 1; idx < argc; idx++)
  {
    inputfiles.emplace_back(argv[idx]);
  }
  if (guess_if_protein_file(inputfiles))
  {
    exit(EXIT_SUCCESS);
  }
  exit(EXIT_FAILURE);
}

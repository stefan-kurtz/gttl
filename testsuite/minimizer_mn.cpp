/*
'\n'ht (c) 2021 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2021 Center for Bioinformatics, University of Hamburg

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
#include <cassert>
#include <cinttypes>
#include <climits>
#include <cstdlib>
#include <cstdint>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <tuple>
#include <vector>

#include "sequences/gttl_multiseq.hpp"
#include "sequences/inputfiles_multiseq.hpp"
#include "sequences/qgrams_hash_nthash.hpp"
#include "sequences/qgrams_hash_nthash.hpp"
#include "sequences/hashed_qgrams.hpp"
#include "minimizer_opt.hpp"
#include "utilities/str_format.hpp"
#include "utilities/runtime_class.hpp"
#include "utilities/constexpr_for.hpp"

std::pair<int,int> determine_hash_bits(int sequences_bits,
                                       int requested_hash_bits)
{
  for (int bytes = 8; bytes <= 9; bytes++)
  {
    const int bits = bytes * CHAR_BIT;
    if (sequences_bits + requested_hash_bits <= bits)
    {
      return std::make_pair(bits - sequences_bits,bytes);
    }
  }
  StrFormat msg("cannot handle sequences_bits + hash_bits = %d + %d > 72",
                sequences_bits,requested_hash_bits);
  throw std::runtime_error{msg.str()};
}

void run_nt_minimizer(const MinimizerOptions &options)
{
  RunTimeClass rt_create_multiseq{};
  constexpr const bool store_header = true;
  constexpr const bool store_sequence = true;
  constexpr const uint8_t padding_char = UINT8_MAX;
  GttlMultiseq *multiseq = gttl_inputfiles_multiseq(options.inputfiles_get(),
                                                    store_header,
                                                    store_sequence,
                                                    padding_char,
                                                    options
                                                      .canonical_option_is_set()
                                                   );
  rt_create_multiseq.show("reading input files and creating multiseq");
  for (auto &log : multiseq->statistics())
  {
    std::cout << "# " << log << '\n';
  }
  int hash_bits = -1;
  int var_sizeof_unit_hashed_qgram;
  assert(options.hash_bits_get() != -1);
  try
  {
    std::tie(hash_bits,var_sizeof_unit_hashed_qgram)
      = determine_hash_bits(multiseq->sequences_bits_get(),
                            options.hash_bits_get());
  }
  catch (const std::exception &err)
  {
    delete multiseq;
    throw;
  }
  assert(var_sizeof_unit_hashed_qgram == 8 or
         var_sizeof_unit_hashed_qgram == 9);
  constexpr_for<8,9+1,1>([&](auto sizeof_unit_hashed_qgram)
  {
    if (sizeof_unit_hashed_qgram == var_sizeof_unit_hashed_qgram)
    {
      std::vector<std::string> log_vector{};
      if (options.canonical_option_is_set())
      {
        using HashedQgrams = HashedQgramsGeneric<sizeof_unit_hashed_qgram,
                                                 QgramNtHashIterator4>;
        HashedQgrams hqg (*multiseq,
                          options.number_of_threads_get(),
                          options.qgram_length_get(),
                          options.window_size_get(),
                          hash_bits,
                          options.sort_by_hash_value_option_is_set(),
                          options.at_constant_distance_option_is_set(),
                          &log_vector);
        RunTimeClass rt_output_hashed_qgrams{};
        hqg.show();
        log_vector.push_back(rt_output_hashed_qgrams
                             .to_string("output of hashed kmers"));
      } else
      {
        using HashedQgrams = HashedQgramsGeneric<sizeof_unit_hashed_qgram,
                                                 QgramNtHashFwdIterator4>;
        HashedQgrams hqg (*multiseq,
                          options.number_of_threads_get(),
                          options.qgram_length_get(),
                          options.window_size_get(),
                          hash_bits,
                          options.sort_by_hash_value_option_is_set(),
                          options.at_constant_distance_option_is_set(),
                          &log_vector);
        if (options.show_mode_get() != 0)
        {
          RunTimeClass rt_output_hashed_qgrams{};
          if (options.show_mode_get() == 1)
          {
            hqg.show();
          } else
          {
            assert(options.show_mode_get() == 2);
            for (auto &&dhqg : hqg)
            {
              printf("%" PRIu64 "\t%zu\t%zu\n",
                     dhqg.hash_value,
                     dhqg.sequence_number,
                     dhqg.startpos);
            }
          }
          log_vector.push_back(rt_output_hashed_qgrams
                               .to_string("output of hashed kmers"));
        }
      }
      for (auto &msg : log_vector)
      {
        std::cout << "# " << msg << '\n';
      }
    }
  });
  delete multiseq;
}

int main(int argc, char *argv[])
{
  RunTimeClass rt_total{};
  MinimizerOptions options{};

  try
  {
    options.parse(argc, argv);
  }
  catch (const std::invalid_argument &err) /* check_err.py */
  {
    std::cerr << argv[0] << ": file " << err.what() << '\n';
    return EXIT_FAILURE;
  }
  if (options.help_option_is_set())
  {
    return EXIT_SUCCESS;
  }
  try
  {
    run_nt_minimizer(options);
  }
  catch (const std::exception &err)
  {
    for (auto &&inputfile : options.inputfiles_get())
    {
      std::cerr << argv[0] << ": file \"" << inputfile << "\"" << err.what()
                << '\n';
    }
    return EXIT_FAILURE;
  }
  rt_total.show("overall");
  return EXIT_SUCCESS;
}

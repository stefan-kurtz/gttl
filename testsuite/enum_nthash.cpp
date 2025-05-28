/*
  Copyright (c) 2021 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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
#include <cstdlib>
#include <cstdint>
#include <exception>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <iostream>
#include <tuple>
#include <cinttypes>
#include "sequences/dna_seq_encoder.hpp"
#include "sequences/gttl_fasta_generator.hpp"
#include "utilities/str_format.hpp"
#include "utilities/mathsupport.hpp"
#include "utilities/cxxopts.hpp"
#include "utilities/unused.hpp"
#include "utilities/runtime_class.hpp"
#include "utilities/bytes_unit.hpp"
#include "sequences/char_range.hpp"
#include "sequences/char_finder.hpp"
#include "sequences/qgrams_hash_nthash.hpp"
#include "sequences/guess_if_protein_seq.hpp"
#include "sequences/gttl_multiseq.hpp"
#include "utilities/wyhash.hpp"

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << '\n';
}

class NtHashOptions
{
 private:
  std::vector<std::string> inputfiles;
  bool help_option,
       bytes_unit_option,
       with_rc_option,
       show_hash_values;
  int hashbits;
  size_t kmer_length;

 public:
  NtHashOptions(void)
    : inputfiles({})
    , help_option(false)
    , bytes_unit_option(false)
    , with_rc_option(false)
    , hashbits(0)
    , kmer_length(0)
  {};

  void parse(int argc, char **argv)
  {
    cxxopts::Options options(argv[0],"process DNA sequence files and "
                                     "generate nthash codes");
    options.set_width(80);
    options.custom_help(std::string("[options] filename0 [filename1]"));
    options.set_tab_expansion();
    options.add_options()
      ("k,kmer_length", "k-mer length",
       cxxopts::value<size_t>(kmer_length)->default_value("18"))
      ("b,hashbits", "number of bits used for hashing",
       cxxopts::value<int>(hashbits)->default_value("39"))
      ("u,bytes_unit", "create BytesUnit object for each kmer and "
                        "corresponding positional information",
        cxxopts::value<bool>(bytes_unit_option)->default_value("false"))
      ("with_rc", "also compute hash values of reverse complement",
       cxxopts::value<bool>(with_rc_option)->default_value("false"))
      ("s,show_hash_values", "show the hash values",
       cxxopts::value<bool>(show_hash_values)->default_value("false"))
     ("h,help", "print usage");
    try
    {
      auto result = options.parse(argc, argv);
      if (result.count("help") > 0)
      {
        help_option = true;
        usage(options);
      } else
      {
        const std::vector<std::string>& unmatched_args = result.unmatched();
        for (size_t idx = 0; idx < unmatched_args.size(); idx++)
        {
          inputfiles.push_back(unmatched_args[idx]);
        }
        if (inputfiles.size() < 1)
        {
          throw cxxopts::exceptions::exception("not enough input files");
        }
      }
    }
    catch (const cxxopts::exceptions::exception &e)
    {
      usage(options);
      throw std::invalid_argument(e.what());
    }
  }
  bool help_option_is_set(void) const noexcept
  {
    return help_option;
  }
  int hashbits_get(void) const noexcept
  {
    return hashbits;
  }
  size_t kmer_length_get(void) const noexcept
  {
    return kmer_length;
  }
  bool bytes_unit_option_is_set(void) const noexcept
  {
    return bytes_unit_option;
  }
  bool with_rc_option_is_set(void) const noexcept
  {
    return with_rc_option;
  }
  bool show_hash_values_is_set(void) const noexcept
  {
    return show_hash_values;
  }
  const std::vector<std::string> &inputfiles_get(void) const noexcept
  {
    return inputfiles;
  }
};

template<class HashValueIterator,
         bool with_rc,
         bool show_hash_values,
         int sizeof_unit_hashed_qgrams,
         bool create_bytes_unit>
static std::tuple<uint64_t,uint64_t,size_t,size_t,size_t> apply_qgram_iterator(
                      size_t qgram_length,
                      uint64_t hashmask,
                      GTTL_UNUSED
                      const GttlBitPacker<sizeof_unit_hashed_qgrams,3>
                        *hashed_qgram_packer,
                      GTTL_UNUSED size_t seqnum,
                      const char *substring,
                      size_t this_length)
{
  HashValueIterator qgiter(qgram_length,substring,this_length);
#ifndef NDEBUG
  uint8_t *qgram_buffer = new uint8_t [qgram_length];
  NThashTransformer nt_hash_transformer(qgram_length);
  auto alphabet = HashValueIterator::alphabet;
#endif

  uint64_t sum_hash_values = 0;
  GTTL_UNUSED uint64_t sum_rc_hash_values = 0;
  size_t seqpos = 0;
  GTTL_UNUSED size_t bytes_unit_sum = 0;
  GTTL_UNUSED size_t bytes_unit_sum_rc = 0;
  for (auto const &&code_pair : qgiter)
  {
    uint64_t this_hash = std::get<0>(code_pair);
    sum_hash_values += this_hash;
#ifndef NDEBUG
    for (size_t idx = 0; idx < qgram_length; idx++)
    {
      qgram_buffer[idx] = alphabet.char_to_rank(substring[seqpos + idx]);
    }
    assert(this_hash == nt_hash_transformer.first_fwd_hash_value_get(
                                                                 qgram_buffer,
                                                                 qgram_length));
#endif
    if constexpr (create_bytes_unit)
    {
      BytesUnit<sizeof_unit_hashed_qgrams,3>
               current_hashed_qgram(*hashed_qgram_packer,
                                    {this_hash & hashmask,
                                     static_cast<uint64_t>(seqnum),
                                     static_cast<uint64_t>(seqpos)});
      bytes_unit_sum += current_hashed_qgram.sum();
    }
    if constexpr (with_rc)
    {
      uint64_t this_rc_hash = std::get<1>(code_pair);
#ifndef NDEBUG
      for (size_t idx = 0; idx < qgram_length; idx++)
      {
        qgram_buffer[qgram_length - 1 - idx]
          = complement_uint8(alphabet.char_to_rank(substring[seqpos + idx]));
      }
      assert(this_rc_hash ==
             nt_hash_transformer.first_fwd_hash_value_get(qgram_buffer,
                                                          qgram_length));
#endif
      if constexpr (show_hash_values)
      {
        printf("%" PRIu64 "\t%" PRIu64 "\n",this_hash,this_rc_hash);
      }
      sum_rc_hash_values += this_rc_hash;
      if constexpr (create_bytes_unit)
      {
        BytesUnit<sizeof_unit_hashed_qgrams,3>
                 current_rc_hashed_qgram(*hashed_qgram_packer,
                                         {this_rc_hash & hashmask,
                                          static_cast<uint64_t>(seqnum),
                                          static_cast<uint64_t>(seqpos)});
        bytes_unit_sum_rc += current_rc_hashed_qgram.sum();
      }
    } else
    {
      if constexpr (show_hash_values)
      {
        printf("%" PRIu64 "\n",this_hash);
      }
    }
    seqpos++;
  }
#ifndef NDEBUG
  delete [] qgram_buffer;
#endif
  return {sum_hash_values,sum_rc_hash_values,seqpos,
          bytes_unit_sum,bytes_unit_sum_rc};
}

template<class HashValueIterator,
         bool with_rc,
         bool show_hash_values,
         int sizeof_unit_hashed_qgrams,
         bool create_bytes_unit,
         bool is_aminoacid>
static void enumerate_nt_hash_template(const char *inputfilename,
                                       size_t qgram_length,
                                       int hashbits,
                                       GTTL_UNUSED int sequences_number_bits,
                                       GTTL_UNUSED int sequences_length_bits)
{
  GttlFpType in_fp = gttl_fp_type_open(inputfilename,"rb");
  if (in_fp == nullptr)
  {
    throw std::runtime_error(": cannot open file");
    /* check_err.py checked */
  }
  size_t total_length = 0;
  size_t seqnum = 0;
  size_t max_sequence_length = 0;
  uint64_t sum_hash_values = 0;
  uint64_t sum_rc_hash_values = 0;
  size_t count_all_qgrams = 0;
  size_t bytes_unit_sum = 0;
  size_t bytes_unit_sum_rc = 0;

  GttlBitPacker<sizeof_unit_hashed_qgrams,3> *hashed_qgram_packer = nullptr;
  if constexpr (create_bytes_unit)
  {
    hashed_qgram_packer = new GttlBitPacker<sizeof_unit_hashed_qgrams,3>
                                           ({hashbits,sequences_number_bits,
                                             sequences_length_bits});
  }
  const uint64_t hashmask = gttl_bits2maxvalue<uint64_t>(hashbits);
  constexpr const int buf_size = 1 << 14;
  GttlFastAGenerator<buf_size> fasta_gen(in_fp);


  for (auto &&si : fasta_gen)
  {
    auto sequence = si->sequence_get();
    total_length += sequence.size();
    max_sequence_length = std::max(max_sequence_length,sequence.size());
    NtCardRanger<is_aminoacid> nuc_ranger(sequence.data(),sequence.size());
    for (auto const &&range : nuc_ranger)
    {
      const size_t this_length = std::get<1>(range);
      const char *substring = sequence.data() + std::get<0>(range);
      auto result = apply_qgram_iterator<HashValueIterator,
                                         with_rc,
                                         show_hash_values,
                                         sizeof_unit_hashed_qgrams,
                                         create_bytes_unit>
                                        (qgram_length,
                                         hashmask,
                                         hashed_qgram_packer,
                                         seqnum,
                                         substring,
                                         this_length);
      sum_hash_values += std::get<0>(result);
      sum_rc_hash_values += std::get<1>(result);
      count_all_qgrams += std::get<2>(result);
      bytes_unit_sum += std::get<3>(result);
      bytes_unit_sum_rc += std::get<4>(result);
    }
    seqnum++;
  }
  delete hashed_qgram_packer;
  printf("# num_of_sequences\t%zu\n",seqnum);
  printf("# total_length\t%zu\n",total_length);
  printf("# max_sequence_length\t%zu\n",max_sequence_length);
  printf("# num_of_sequences.bits\t%d\n",gttl_required_bits(seqnum));
  printf("# max_sequence_length.bits\t%d\n",
          gttl_required_bits<size_t>(max_sequence_length));
  printf("# count_all_qgrams\t%zu\n",count_all_qgrams);
  printf("# sum_hash_values\t%" PRIu64 "\n",sum_hash_values);
  printf("# bytes_unit_sum\t%zu\n", bytes_unit_sum);
  if constexpr (with_rc)
  {
    printf("# sum_rc_hash_values\t%" PRIu64 "\n",sum_rc_hash_values);
    printf("# bytes_unit_sum_rc\t%zu\n",bytes_unit_sum_rc);
  }
}

#define CALL_enumerate_nt_hash_template(BYTE_UNITS,BYTE_UNITS_OPTION)\
        if (show_hash_values)\
        {\
          enumerate_nt_hash_template<HashValueIterator,\
                                     with_rc,\
                                     true,\
                                     BYTE_UNITS,BYTE_UNITS_OPTION,\
                                     is_aminoacid>\
                                    (inputfilename,\
                                     qgram_length,\
                                     hashbits,\
                                     sequences_number_bits,\
                                     sequences_length_bits);\
        } else\
        {\
          enumerate_nt_hash_template<HashValueIterator,\
                                     with_rc,\
                                     false,\
                                     BYTE_UNITS,BYTE_UNITS_OPTION,\
                                     is_aminoacid>\
                                    (inputfilename,\
                                     qgram_length,\
                                     hashbits,\
                                     sequences_number_bits,\
                                     sequences_length_bits);\
        }

template<class HashValueIterator,bool with_rc, bool is_aminoacid>
static void enumerate_nt_hash(const char *inputfilename,
                              bool show_hash_values,
                              bool bytes_unit_option,
                              size_t qgram_length,
                              int hashbits)
{
  constexpr const bool store_header = false;
  constexpr const bool store_sequence = false;
  GttlMultiseq multiseq(inputfilename,store_header,store_sequence,UINT8_MAX);
  const int sequences_number_bits = multiseq.sequences_number_bits_get();
  const int sequences_length_bits = multiseq.sequences_length_bits_get();
  RunTimeClass rt_nthash{};
  if (sequences_number_bits + sequences_length_bits + hashbits <= 64)
  {
    if (bytes_unit_option)
    {
      CALL_enumerate_nt_hash_template(8,true);
    } else
    {
      CALL_enumerate_nt_hash_template(8,false);
    }
  } else
  {
    if (bytes_unit_option)
    {
      CALL_enumerate_nt_hash_template(9,true);
    } else
    {
      CALL_enumerate_nt_hash_template(9,false);
    }
  }
  StrFormat msg("nthash\t%s\t%zu\t%d\t%d\t%d",
                inputfilename,qgram_length,
                hashbits,sequences_number_bits,
                sequences_length_bits);
  rt_nthash.show(msg.str());
}

int main(int argc,char *argv[])
{
  NtHashOptions options{};

  try
  {
    options.parse(argc, argv);
  }
  catch (std::invalid_argument &e) /* check_err.py */
  {
    std::cerr << argv[0] << ": " << e.what() << '\n';
    return EXIT_FAILURE;
  }
  if (options.help_option_is_set())
  {
    return EXIT_SUCCESS;
  }
  const char *progname = argv[0];
  bool haserr = false;
  for (auto &&inputfile : options.inputfiles_get())
  {
    try
    {
      const bool is_protein = guess_if_protein_file(inputfile.c_str());
      if (is_protein)
      {
        if (options.with_rc_option_is_set())
        {
          // There is no such thing as a reverse-complement
          // of an aminoacid sequence
          std::cerr <<
            "Cannot run reverse-complements on aminoacid sequences!\n";
          exit(EXIT_FAILURE);
        }

        RunTimeClass rt{};
        enumerate_nt_hash<QgramNtHashAAFwdIterator20, false, true>
                         (inputfile.c_str(),
                          options.show_hash_values_is_set(),
                          options.bytes_unit_option_is_set(),
                          options.kmer_length_get(),
                          options.hashbits_get());
        rt.show("NtHash with input: ");

#ifdef WYHASH_STATISTICS
        GttlFpType fp = gttl_fp_type_open(inputfile.c_str(), "rb");

        constexpr const int buf_size = 1 << 14;
        GttlFastAGenerator<buf_size> fasta_gen(fp);


        using NucleotideRanger
          = GttlCharRange<typename RangerTraits<true>::Finder,
                          RangerTraits<true>::instance,
                          true,false>;

        rt.reset();
        size_t total_length = 0;
        size_t seqnum = 0;
        size_t max_sequence_length = 0;
        uint64_t sum_hash_values = 0;
        size_t count_all_qgrams = 0;

        for (auto &&si : fasta_gen)
        {
          auto sequence = si->sequence_get();
          total_length += sequence.size();
          max_sequence_length = std::max(max_sequence_length,sequence.size());
          NucleotideRanger nuc_ranger(sequence.data(),sequence.size());
          for (auto const &&range : nuc_ranger)
          {
            const size_t this_length = std::get<1>(range);
            const char *substring = sequence.data() + std::get<0>(range);

            if(this_length >= options.kmer_length_get())
            {
              for(size_t i = 0;
                  i < this_length - options.kmer_length_get() + 1; i++)
              {
                if(std::strlen(substring + i) >= options.kmer_length_get())
                {
                  const uint64_t hash = wyhash(substring + i,
                                               options.kmer_length_get(),
                                               0xABCDEFABCDEF);
                  sum_hash_values += hash;
                  count_all_qgrams++;
                }
              }
            }
          }
          seqnum++;
        }
        printf("# num_of_sequences\t%zu\n",seqnum);
        printf("# total_length\t%zu\n",total_length);
        printf("# max_sequence_length\t%zu\n",max_sequence_length);
        printf("# num_of_sequences.bits\t%d\n",gttl_required_bits(seqnum));
        printf("# max_sequence_length.bits\t%d\n",
               gttl_required_bits<size_t>(max_sequence_length));
        printf("# count_all_qgrams\t%zu\n",count_all_qgrams);
        printf("# sum_hash_values\t%" PRIu64 "\n",sum_hash_values);
        rt.show("Runtime wyhash: ");
#endif
      } else
      {
        if (options.with_rc_option_is_set())
        {
          enumerate_nt_hash<QgramNtHashIterator4,true, false>
                           (inputfile.c_str(),
                            options.show_hash_values_is_set(),
                            options.bytes_unit_option_is_set(),
                            options.kmer_length_get(),
                            options.hashbits_get());
        } else
        {
          enumerate_nt_hash<QgramNtHashFwdIterator4,false, false>
                           (inputfile.c_str(),
                            options.show_hash_values_is_set(),
                            options.bytes_unit_option_is_set(),
                            options.kmer_length_get(),
                            options.hashbits_get());
        }
      }
    }
    catch (std::exception &msg)
    {
      std::cerr << progname << ": file \"" << inputfile << "\""
                << msg.what() << '\n';
      haserr = true;
      break;
    }
  }
  return haserr ? EXIT_FAILURE : EXIT_SUCCESS;
}

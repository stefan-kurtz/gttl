/*
  Copyright (c) 2013-2022 Stefan Kurtz <stefan.kurtz@uni-hamburg.de>
  Copyright (c) 2013-2022 Center for Bioinformatics, University of Hamburg

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
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <exception>
#include <ios>
#include <string>
#include <ostream>
#include <iostream>
#include <stdexcept>
#include <filesystem>
#include <vector>

#include "sequences/alphabet.hpp"
#include "utilities/bitpacker.hpp"
#include "utilities/bytes_unit.hpp"
#include "utilities/constexpr_for.hpp"
#include "utilities/gttl_binary_read.hpp"
#include "utilities/mathsupport.hpp"
#include "utilities/ordered_integer_sequence.hpp"
#include "utilities/runtime_class.hpp"
#include "utilities/gttl_file_open.hpp"
#include "utilities/gttl_binary_write.hpp"
#include "utilities/memory_tracker.hpp"
#include "sequences/guess_if_protein_seq.hpp"
#include "sequences/gttl_multiseq.hpp"
#include "sequences/literate_multiseq.hpp"
#include "sequences/inputfiles_multiseq.hpp"
#include "indexes/succinct_bitvector.hpp"

#include "sa_induced_options.hpp"
#include "fill_bu_suftab.hpp"
#include "kasai.hpp"
#include "sk_sain.hpp"
#include "plcp_table.hpp"

static void prj_file_write(const std::string &filename,
                           bool reverse_complement,
                           size_t nonspecial_suffixes,
                           size_t sequences_number,
                           int sequences_number_bits,
                           int sequences_length_bits,
                           size_t sizeof_suftab_entry,
                           const std::vector<std::string> &inputfiles)
{
  std::ofstream out_stream;
  out_stream.open(filename);
  out_stream << "reverse_complement\t" << (reverse_complement ? 1 : 0) << '\n';
  out_stream << "nonspecial_suffixes\t" << nonspecial_suffixes << '\n';
  out_stream << "sequences_number\t" << sequences_number << '\n';
  out_stream << "sequences_number_bits\t" << sequences_number_bits << '\n';
  out_stream << "sequences_length_bits\t" << sequences_length_bits << '\n';
  out_stream << "sizeof_suftab_entry\t" << sizeof_suftab_entry << '\n';
  for (auto &&inputfile : inputfiles)
  {
    out_stream << "inputfile\t" << inputfile << '\n';
  }
  out_stream.close();
}

static FILE *outstream_create(const std::string &indexname,
                              const std::string &filenamesuffix)
{
  const std::string outfilename{indexname + filenamesuffix};

  FILE *out_fp = fopen(outfilename.c_str(), "wb");
  if (out_fp == nullptr)
  {
    throw std::ios_base::failure(std::string(": cannot create file \"") +
                                 outfilename + std::string("\""));
  }
  return out_fp;
}

template <typename SuftabBaseType>
static void output_suftab(const SuftabBaseType *suftab,
                          size_t totallength,
                          const std::string &indexname,
                          const std::string &filenamesuffix)
{
  FILE *out_fp = outstream_create(indexname, filenamesuffix);
  std::fwrite(suftab, sizeof *suftab, totallength + 1, out_fp);
  fclose(out_fp);
}

template<class Generator>
static void lcptab_output_saturated(const std::string &indexname,
                                    Generator &lcp_generator)
{
  BinaryFileWriter<uint8_t> lcp_writer(indexname + ".lcp");
  BinaryFileWriter<uint16_t> ll2_writer(indexname + ".ll2");
  BinaryFileWriter<uint32_t> ll4_writer(indexname + ".ll4");
  for (auto &&lcp_value : lcp_generator)
  {
    if (lcp_value < UINT8_MAX)
    {
      lcp_writer.append(static_cast<uint8_t>(lcp_value));
    } else
    {
      lcp_writer.append(static_cast<uint8_t>(UINT8_MAX));
      if (lcp_value < UINT16_MAX)
      {
        ll2_writer.append(static_cast<uint16_t>(lcp_value));
      } else
      {
        ll2_writer.append(static_cast<uint16_t>(UINT16_MAX));
        ll4_writer.append(static_cast<uint32_t>(lcp_value));
      }
    }
  }
}

template<class SuftabBaseType>
static void lcptab_output_succinct(const std::string &indexname,
                                   PlcpTable<SuftabBaseType>
                                     &plcp_table)
{
  SuccinctBitvector b;
  size_t last = size_t(1);
  for (size_t pos = 0; pos < plcp_table.get_total_length(); pos++){
    const size_t cur = plcp_table.plcp_value_get(pos);
    const size_t unary = cur - last + 1;
    /* Implement method to construct unary coding by single function call */
    for (size_t i = 0; i < unary; i++)
    {
      b.push(false);
    }
    // printf("%zu; %zu %zu %zu %zu\n", pos, cur, cur + pos, cur - last + 1,
    //         b.get_length());

    b.push(true);
    last = cur;
  }

  b.buildAccelerationStructures();
  // printf("%zu, %zu\n", plcp_table.get_total_length(),
  //                      b.get_rank(b.get_length(), 1));
  const std::string lls_filename(indexname + ".lls");
  b.serialize(lls_filename);
  // b.print();
}

template<class Generator>
static void lcptab_print(Generator &lcp_generator)
{
  std::cout << "# lcptab: \n";
  for(auto &&lcp_value : lcp_generator)
  {
    std::cout << lcp_value << std::endl;
  }
}

template<typename SuftabBaseType>
static void create_and_output_lcptab(GttlMemoryTracker *memory_tracker,
                                     const SuftabBaseType *suftab,
                                     const std::string &indexname,
                                     const uint8_t *seq,
                                     size_t totallength,
                                     size_t var_alphasize,
                                     bool print_stdout,
                                     SainOptions::LcptabMethod lcptab_method,
                                     bool succinct_lcptab)
{
  if (suftab != nullptr or
      lcptab_method == SainOptions::LcptabMethod::Lcptab_kasai9n)
  {
    std::vector<uint32_t> lcptab
      = suftab != nullptr
          ? gttl_lcp13_kasai<SuftabBaseType>(memory_tracker,
                                             seq,
                                             totallength,
                                             suftab,
                                             var_alphasize)
          : gttl_lcp9_kasai<SuftabBaseType>(memory_tracker,
                                            seq,
                                            totallength,
                                            indexname,
                                            var_alphasize);
    assert(lcptab[0] == 0);
    assert(lcptab[totallength] == 0);
    lcptab_output_saturated(indexname, lcptab);
    if (print_stdout)
    {
      lcptab_print(lcptab);
    }
    memory_tracker->untrack(lcptab.data(),__FILE__,__LINE__);
  } else
  {
    assert (lcptab_method == SainOptions::LcptabMethod::Lcptab_plcp5n);
    PlcpTable<SuftabBaseType> plcp_table(memory_tracker,
                                         seq, totallength,
                                         indexname, var_alphasize);
    if (succinct_lcptab)
    {
      lcptab_output_succinct(indexname, plcp_table);
    } else
    {
      lcptab_output_saturated(indexname, plcp_table);
    }
    if (print_stdout)
    {
      lcptab_print(plcp_table);
    }
  }
}

template<typename SuftabBaseType>
static void show_suftab(const SuftabBaseType *suftab, size_t totallength)
{
  for (size_t idx = 0; idx <= totallength; idx++)
  {
    std::cout << suftab[idx] << std::endl;
  }
}

template <typename T_suf, typename T_bp>
static void show_bu_suftab(const std::string &indexname, const T_bp &bp)
{
  std::cout << "# Suftab:\n";
  const std::string bu_suftab_filename{indexname + ".bsf"};
  const BinaryFileReader<T_suf> bu_suftab_reader(bu_suftab_filename);
  for (auto &bu_suftab_entry : bu_suftab_reader)
  {
    std::cout << bu_suftab_entry.template decode_at<0>(bp) << "\t";
    std::cout << bu_suftab_entry.template decode_at<1>(bp) << std::endl;
  }
}

template<typename SuftabBaseType>
static void output_tis_suf_tables(const SainOptions &sainoptions,
                                  const SuftabBaseType *suftab,
                                  const uint8_t *filecontents,
                                  size_t totallength)
{
  const std::string &indexname = sainoptions.indexname_get();
  if (sainoptions.tistab_show_opt_is_set())
  {
    std::cout << "# tistab: \n";
    std::fwrite(filecontents, sizeof *filecontents, totallength, stdout);
  }
  if (sainoptions.tistab_out_opt_is_set())
  {
    FILE *out_fp = outstream_create(indexname, ".tis");
    std::fwrite(filecontents, sizeof *filecontents, totallength, out_fp);
    fclose(out_fp);
  }
  if (sainoptions.abs_suftab_show_opt_is_set())
  {
    show_suftab<SuftabBaseType>(suftab, totallength);
  }
  if (sainoptions.abs_suftab_out_opt_is_set())
  {
    output_suftab<SuftabBaseType>(suftab, totallength, indexname, ".suf");
  }
}

template<typename SuftabBaseType>
static void enhanced_suffixarray_plain_input_format(GttlMemoryTracker
                                                      *memory_tracker,
                                                    const SainOptions
                                                      &sainoptions,
                                                    const uint8_t *filecontents,
                                                    size_t totallength)
{
  const std::string &indexname = sainoptions.indexname_get();
  constexpr const bool intermediatecheck = false;
  const int required_bits = gttl_required_bits(totallength);
  static constexpr const size_t const_alphasize = UINT8_MAX + 1;
  SuftabBaseType *suftab
    = gttl_sain_plain_sortsuffixes<SuftabBaseType,const_alphasize>
                                  (memory_tracker,
                                   sainoptions.verbose_opt_is_set(),
                                   filecontents,
                                   totallength,
                                   intermediatecheck,
                                   sainoptions.check_suftab_opt_is_set());
  try
  {
    output_tis_suf_tables<SuftabBaseType>(sainoptions,
                                          suftab,
                                          filecontents,
                                          totallength);
  }
  catch (const std::exception &err)
  {
    GTTL_UNTRACK_ALLOC(suftab);
    free(suftab);
    suftab = nullptr;
    throw;
  }
  RunTimeClass rt_construct_lcp{};
  if (sainoptions.lcptab_method_get() != SainOptions::Lcptab_no)
  {
    create_and_output_lcptab<SuftabBaseType>(memory_tracker,
                                             suftab,
                                             indexname,
                                             filecontents,
                                             totallength,
                                             const_alphasize,
                                             sainoptions
                                               .lcptab_show_opt_is_set(),
                                             sainoptions.lcptab_method_get(),
                                             sainoptions
                                               .succinct_option_is_set());
  }
  if (sainoptions.verbose_opt_is_set() and
      (sainoptions.lcptab_method_get() != SainOptions::Lcptab_no or
       sainoptions.lcptab_show_opt_is_set()))
  {
    rt_construct_lcp.show("construction of lcp table");
  }
  GTTL_UNTRACK_ALLOC(suftab);
  free(suftab);
  suftab = nullptr;
  const size_t nonspecial_suffixes = totallength;
  const size_t sequences_number = 1;
  const int sequences_number_bits = 0;
  const int sequences_length_bits = required_bits;
  prj_file_write(indexname + ".prj",
                 sainoptions.reverse_complement_option_is_set(),
                 nonspecial_suffixes,
                 sequences_number,
                 sequences_number_bits,
                 sequences_length_bits,
                 sizeof(SuftabBaseType) * CHAR_BIT,
                 sainoptions.inputfiles_get());
}

template<typename SuftabBaseType>
static void enhanced_suffixarray_multiseq(GttlMemoryTracker *memory_tracker,
                                          const SainOptions &sainoptions,
                                          GttlMultiseq *multiseq,
                                          size_t totallength,
                                          bool is_protein)
{
  const std::string &indexname = sainoptions.indexname_get();
  constexpr const bool intermediatecheck = false;
  const int sequences_number_bits = multiseq->sequences_number_bits_get();
  const int sequences_length_bits = multiseq->sequences_length_bits_get();
  const int required_bits = sequences_number_bits + sequences_length_bits;
  memory_tracker->track(multiseq,__FILE__,__LINE__,multiseq->size_in_bytes());

  if (sainoptions.verbose_opt_is_set())
  {
    std::cout << "# totallength\t" << totallength << '\n';
    std::cout << "# number of sequences\t" << multiseq->sequences_number_get()
              << '\n';
    std::cout << "# alphabet size\t" << (is_protein ? 20 : 4) << '\n';
    std::cout << "# sequences_number_bits\t" << sequences_number_bits << '\n';
    std::cout << "# sequences_length_bits\t" << sequences_length_bits << '\n';
    std::cout << "# required_bits\t" << required_bits << '\n';
    std::cout << "# SPACE\tsequences representation (MB):\t"
              << static_cast<size_t>(mega_bytes(multiseq->size_in_bytes()))
              << '\n';
  }

  SuftabBaseType *suftab = nullptr;
  size_t nonspecial_suffixes = 0;
  constexpr_for<4, 20 + 1, 16>([&](auto const_alphasize) {
    if (const_alphasize == 4 and not is_protein)
    {
      LiterateMultiseq<alphabet::nucleotides_upper_lower,const_alphasize>
         lit_multiseq (multiseq);
      lit_multiseq.perform_sequence_encoding();
      size_t *charcount_arr = lit_multiseq.rank_dist_copy();
      memory_tracker->track(charcount_arr,__FILE__,__LINE__,
                            lit_multiseq.number_of_ranks() *
                            sizeof *charcount_arr);
      const size_t numofwildcards = charcount_arr[const_alphasize];
      nonspecial_suffixes
        = multiseq->sequences_total_length_get() - numofwildcards;
      if (sainoptions.verbose_opt_is_set())
      {
        std::cout << "# special suffixes\t" << numofwildcards << '\n';
        std::cout << "# nonspecial suffixes\t" << nonspecial_suffixes << '\n';
      }
      suftab = gttl_sain_multiseq_sortsuffixes<SuftabBaseType,const_alphasize>
                                              (memory_tracker,
                                               sainoptions.verbose_opt_is_set(),
                                               multiseq,
                                               charcount_arr,
                                               intermediatecheck,
                                               sainoptions
                                                 .check_suftab_opt_is_set(),
                                               sainoptions
                                                 .buffered_option_is_set());
      memory_tracker->untrack(charcount_arr,__FILE__,__LINE__);
      free(charcount_arr);
    } else
    {
      if (const_alphasize == 20 and is_protein)
      {
        LiterateMultiseq<alphabet::amino_acids, const_alphasize>
                      lit_multiseq(multiseq);
        lit_multiseq.perform_sequence_encoding();
        size_t *charcount_arr = lit_multiseq.rank_dist_copy();
        memory_tracker->track(charcount_arr,__FILE__,__LINE__,
                              lit_multiseq.number_of_ranks() *
                              sizeof *charcount_arr);
        const size_t numofwildcards = charcount_arr[const_alphasize];
        nonspecial_suffixes
          = multiseq->sequences_total_length_get() - numofwildcards;
        if (sainoptions.verbose_opt_is_set())
        {
          std::cout << "# special suffixes\t" << numofwildcards << '\n';
          std::cout << "# nonspecial suffixes\t" << nonspecial_suffixes << '\n';
        }
        suftab
          = gttl_sain_multiseq_sortsuffixes<SuftabBaseType,const_alphasize>
                                           (memory_tracker,
                                            sainoptions.verbose_opt_is_set(),
                                            multiseq,
                                            charcount_arr,
                                            intermediatecheck,
                                            sainoptions
                                              .check_suftab_opt_is_set(),
                                            sainoptions
                                              .buffered_option_is_set());
        memory_tracker->untrack(charcount_arr,__FILE__,__LINE__);
        free(charcount_arr);
      }
    }
  });
  const size_t sequences_number = multiseq->sequences_number_get();
  try
  {
    output_tis_suf_tables<SuftabBaseType>
                         (sainoptions,
                          suftab,
                          reinterpret_cast<const uint8_t *>
                            (multiseq->sequence_ptr_get()),
                          totallength);
  }
  catch (const std::exception &err)
  {
    GTTL_UNTRACK_ALLOC(suftab);
    free(suftab);
    suftab = nullptr;
    throw;
  }
  RunTimeClass rt_construct_lcp{};
  if (sainoptions.lcptab_method_get() != SainOptions::LcptabMethod::Lcptab_no)
  {
    if (sainoptions.lcptab_method_get() == SainOptions::Lcptab_kasai9n or
        sainoptions.lcptab_method_get() == SainOptions::Lcptab_plcp5n)
    {
      GTTL_UNTRACK_ALLOC(suftab);
      free(suftab);
      suftab = nullptr;
    }
    create_and_output_lcptab<SuftabBaseType>(memory_tracker,
                                             suftab,
                                             indexname,
                                             reinterpret_cast<const uint8_t *>
                                               (multiseq->sequence_ptr_get()),
                                             totallength,
                                             is_protein ? size_t(20)
                                                        : size_t(4),
                                             sainoptions
                                               .lcptab_show_opt_is_set(),
                                             sainoptions.lcptab_method_get(),
                                             sainoptions
                                               .succinct_option_is_set());
  }
  if (sainoptions.verbose_opt_is_set() and
      (sainoptions.lcptab_method_get() != SainOptions::Lcptab_no or
       sainoptions.lcptab_show_opt_is_set()))
  {
    rt_construct_lcp.show("construction of lcp table");
  }
  if (sainoptions.rel_suftab_out_opt_is_set() or
      sainoptions.rel_suftab_show_opt_is_set())
  {
    constexpr_for<4, 7 + 1, 1>([&](auto req_size)
    {
      if ((req_size == 4 and required_bits <= 32) or
          (required_bits > (req_size - 1) * 8 and
           required_bits <= 8 * req_size))
      {
        static constexpr const int bit_groups = 2;
        static constexpr const int sizeof_unit = req_size;
        using SainBitPacker = GttlBitPacker<sizeof_unit, bit_groups>;
        using SainBytesUnit = BytesUnit<sizeof_unit, bit_groups>;
        int first_group_bits;

        if (multiseq->sequences_number_get() == 1)
        {
          assert(sequences_number_bits == 0);
          first_group_bits = sizeof_unit * CHAR_BIT - sequences_length_bits;
        } else
        {
          assert(sequences_number_bits > 0);
          first_group_bits = sequences_number_bits;
        }
        const SainBitPacker bp({first_group_bits, sequences_length_bits});

        RunTimeClass rt_construct_bu_suftab{};
        if (sequences_number_bits == 0 and req_size == sizeof(SuftabBaseType))
        {
          if (sainoptions.rel_suftab_out_opt_is_set())
          {
            if (sainoptions.abs_suftab_out_opt_is_set())
            {
              const std::string suftab_path(indexname + ".suf");
              const std::string bu_suftab_path(indexname + ".bsf");
              assert(std::filesystem::exists(suftab_path));
              if (std::filesystem::is_symlink(bu_suftab_path) or
                  std::filesystem::exists(bu_suftab_path))
              {
                std::remove(bu_suftab_path.c_str());
              }
              auto suftab_filename
                = std::filesystem::path(suftab_path).filename();
              assert(not suftab_filename.empty());
              std::filesystem::create_symlink(suftab_filename, bu_suftab_path);
            } else
            {
              try
              {
                assert(suftab != nullptr);
                SainBytesUnit *bu_suftab
                  = reinterpret_cast<SainBytesUnit *>(suftab);
                output_suftab<SainBytesUnit>(bu_suftab, totallength, indexname,
                                             ".bsf");
              }
              catch (const std::exception &err)
              {
                GTTL_UNTRACK_ALLOC(suftab);
                free(suftab);
                suftab = nullptr;
                throw;
              }
            }
          }
          if (sainoptions.rel_suftab_show_opt_is_set())
          {
            show_bu_suftab<SainBytesUnit, SainBitPacker>(indexname, bp);
          }
        } else
        {
          if (sainoptions.intset_sizeof_get() == -1)
          {
            bool bu_suftab_own;
            SainBytesUnit *bu_suftab;
            if (sainoptions.abs_suftab_out_opt_is_set())
            {
              bu_suftab_own = true;
              if (not sainoptions.check_suftab_opt_is_set() and
                  suftab != nullptr)
              {
                GTTL_UNTRACK_ALLOC(suftab);
                free(suftab);
                suftab = nullptr;
              }
              bu_suftab = fill_bu_suftab_MULTISEQ_linear<SuftabBaseType,
                                                         SainBitPacker,
                                                         SainBytesUnit>
                                                         (memory_tracker,
                                                          indexname,
                                                          bp,
                                                          multiseq,
                                                          totallength);
            } else
            {
              bu_suftab_own = sainoptions.check_suftab_opt_is_set() or
                              (sizeof *bu_suftab != sizeof *suftab);
              if (bu_suftab_own)
              {
                bu_suftab = new SainBytesUnit [totallength + 1];
                memory_tracker->track(bu_suftab,__FILE__,__LINE__,
                                      (totallength + 1)
                                         * sizeof(SainBytesUnit));
              } else
              {
                assert(suftab != nullptr);
                bu_suftab = reinterpret_cast<SainBytesUnit *>(suftab);
              }
              fill_bu_suftab_MULTISEQ_linear<SuftabBaseType,
                                             SainBitPacker,
                                             SainBytesUnit>
                                             (memory_tracker,
                                              suftab,
                                              bu_suftab,
                                              bp,
                                              multiseq,
                                              totallength);
            }
            FILE *out_fp;
            try
            {
              out_fp = outstream_create(indexname,".bsf");
            }
            catch (const std::exception &err)
            {
              if (bu_suftab_own)
              {
                memory_tracker->untrack(bu_suftab,__FILE__,__LINE__);
                delete[] bu_suftab;
              }
              if (suftab != nullptr)
              {
                GTTL_UNTRACK_ALLOC(suftab);
                free(suftab);
                suftab = nullptr;
              }
              throw;
            }
            std::fwrite(bu_suftab,sizeof *bu_suftab,totallength + 1,out_fp);
            fclose(out_fp);
            if (sainoptions.rel_suftab_show_opt_is_set())
            {
              show_bu_suftab<SainBytesUnit, SainBitPacker>(indexname, bp);
            }
            if (bu_suftab_own)
            {
              memory_tracker->untrack(bu_suftab,__FILE__,__LINE__);
              delete[] bu_suftab;
            }
          } else
          {
            size_t sizeof_intset_choice;
            if (sainoptions.intset_sizeof_get() == 0)
            {
              const size_t sequences_number
                = multiseq->sequences_number_get();
              const size_t totallength
                = multiseq->sequences_total_length_get() + sequences_number
                                                         - 1;
              sizeof_intset_choice
                = ordered_integer_sequence_sizeof_smallest(totallength,
                                                           sequences_number);
            } else
            {
              sizeof_intset_choice = sainoptions.intset_sizeof_get();
            }
            if (sizeof_intset_choice == 1)
            {
#define FILL_BU_SUFTAB_MULTISEQ(BASE_TYPE)\
      if (sainoptions.rel_suftab_out_opt_is_set())\
      {\
        fill_bu_suftab_MULTISEQ<SuftabBaseType,SainBitPacker,SainBytesUnit,\
                                BASE_TYPE>\
                               (suftab,\
                                bp,\
                                multiseq,\
                                totallength,\
                                sainoptions.verbose_opt_is_set(),\
                                indexname);\
      }
              FILL_BU_SUFTAB_MULTISEQ(uint8_t);
            } else
            {
              if (sizeof_intset_choice == 2)
              {
                FILL_BU_SUFTAB_MULTISEQ(uint16_t);
              } else
              {
                assert (sizeof_intset_choice == 4);
                FILL_BU_SUFTAB_MULTISEQ(uint32_t);
              }
            }
          }
        }
        if (sainoptions.verbose_opt_is_set())
        {
          rt_construct_bu_suftab.show("construction of bu_suftab");
        }
        if (sainoptions.check_suftab_opt_is_set())
        {
          check_bu_suftab_MULTISEQ<SuftabBaseType,SainBytesUnit,SainBitPacker>
                                  (suftab, multiseq, bp, totallength,
                                   indexname + ".bsf");
        }
      }
    });
  }
  if (suftab != nullptr)
  {
    GTTL_UNTRACK_ALLOC(suftab);
    free(suftab);
    suftab = nullptr;
  }
  prj_file_write(indexname + ".prj",
                 sainoptions.reverse_complement_option_is_set(),
                 nonspecial_suffixes,
                 sequences_number,
                 sequences_number_bits,
                 sequences_length_bits,
                 sizeof(SuftabBaseType) * CHAR_BIT,
                 sainoptions.inputfiles_get());
}

int main(int argc, char *argv[])
{
  RunTimeClass rt_overall{};
  SainOptions sainoptions{};
  GttlMemoryTracker memory_tracker{};

  try
  {
    sainoptions.parse(argc, argv);
  }
  catch (const std::invalid_argument &err)
  {
    std::cerr << argv[0] << ": " << err.what() << '\n';
    return EXIT_FAILURE;
  }
  if (sainoptions.help_opt_is_set())
  {
    printf("# CXX_VERSION\t%s\n",CXX_VERSION);
    return EXIT_SUCCESS;
  }
  try
  {
    if (sainoptions.verbose_opt_is_set())
    {
      printf("# CXX_VERSION\t%s\n",CXX_VERSION);
      printf("# CXX_FLAGS\t%s\n",CXX_FLAGS);
    }
    constexpr const size_t max_sequence_length_uint32_t = size_t(1) << 30;
    size_t totallength;
    if (sainoptions.plain_input_format_opt_is_set())
    {
      RunTimeClass rt_input_files{};
      auto filecontents
        = gttl_read_files<uint8_t>(sainoptions.inputfiles_get());
      totallength = filecontents.size();
      if (sainoptions.verbose_opt_is_set())
      {
        rt_input_files.show("input of file(s)");
      }
      if (totallength == 0)
      {
        throw std::runtime_error(std::string(": cannot construct suffix array "
                                             "from empty string"));
      }
      memory_tracker.track(filecontents.data(),__FILE__,__LINE__,
                           totallength * sizeof(uint8_t));
      try
      {
        (totallength <= max_sequence_length_uint32_t
          ? enhanced_suffixarray_plain_input_format<uint32_t>
          : enhanced_suffixarray_plain_input_format<uint64_t>)
              (&memory_tracker,
               sainoptions,
               filecontents.data(),
               totallength);
      }
      catch (const std::exception &err)
      {
        memory_tracker.untrack(filecontents.data(),__FILE__,__LINE__);
        throw;
      }
      memory_tracker.untrack(filecontents.data(),__FILE__,__LINE__);
    } else
    {
      GttlMultiseq *multiseq;
      RunTimeClass rt_input_files{};
      const bool is_protein
        = guess_if_protein_file(sainoptions.inputfiles_get());
      if (is_protein and sainoptions.reverse_complement_option_is_set())
      {
        throw std::invalid_argument(
          std::string(": option --reverse_complement is "
                      "only possible for DNA sequences"));
      }
      constexpr const bool store_header = false;
      constexpr const bool store_sequence = true;
      const uint8_t padding_char = is_protein ? uint8_t(20) : uint8_t(4);
      multiseq = gttl_inputfiles_multiseq(sainoptions.inputfiles_get(),
                                          store_header,
                                          store_sequence,
                                          padding_char,
                                          sainoptions
                                            .reverse_complement_option_is_set()
                                         );
      assert(multiseq != nullptr);
      totallength = multiseq->sequences_number_get() +
                    multiseq->sequences_total_length_get() - 1;
      if (sainoptions.verbose_opt_is_set())
      {
        rt_input_files.show("input of file(s)");
      }
      try
      {
        printf("# use uint%d_t as base type for suffix array\n",
                totallength <= max_sequence_length_uint32_t ? 32 : 64);
        (totallength <= max_sequence_length_uint32_t
          ? enhanced_suffixarray_multiseq<uint32_t>
          : enhanced_suffixarray_multiseq<uint64_t>)(&memory_tracker,
                                                     sainoptions,
                                                     multiseq,
                                                     totallength,
                                                     is_protein);
      }
      catch (const std::exception &err)
      {
        memory_tracker.untrack(multiseq,__FILE__,__LINE__);
        delete multiseq;
        throw;
      }
      memory_tracker.untrack(multiseq,__FILE__,__LINE__);
      delete multiseq;
    }
    if (sainoptions.verbose_opt_is_set())
    {
      rt_overall.show("overall");
      printf("# SPACE\tpeak (MB):\t%zu\n",memory_tracker.peak_get());
      const size_t peak_in_bytes = memory_tracker.peak_in_bytes_get();
      printf("# SPACE\tpeak/input symbol (bytes):\t%.2f\n",
              static_cast<double>(peak_in_bytes)/totallength);
    }
  }
  catch (const std::exception &err)
  {
    std::cerr << argv[0] << ": file \"" << sainoptions.inputfiles_get()[0]
              << "\"" << err.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

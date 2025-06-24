#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>
#include "sequences/guess_if_protein_seq.hpp"
#include "sequences/qgrams_hash_nthash.hpp"
#include "utilities/has_fasta_or_fastq_extension.hpp"
#include "utilities/runtime_class.hpp"
#include "utilities/str_format.hpp"
#include "utilities/nttable.hpp"
#include "utilities/binary_nttable.hpp"
#include "sequences/ntcard.hpp"
#include "ntcard_opt.hpp"

template<bool split_at_wildcard>
static void estimate_F_values(const NtcardOptions &options)
{
  static constexpr const uint8_t undefined_rank = 0;
  const bool is_protein =
    gttl_likely_fasta_format(options.inputfile_get())
      ? guess_if_protein_file(options.inputfile_get().c_str())
      : false;

  if (options.binary_option_is_set())
  {
    BinaryNtTable table = is_protein ?
          ntcard_enumerate<split_at_wildcard,
                           QgramNtHashAAFwdIteratorGeneric<undefined_rank>,
                           QgramNtHashAAFwdIteratorGenericNoTransform
                             <undefined_rank>,
                           BinaryNtTable,
                           true>
                          (options.inputfile_get(),
                           options.qgram_length_get(),
                           options.s_get(),
                           options.r_get(),
                           options.num_threads_get())
          :
          ntcard_enumerate<split_at_wildcard,
                           QgramNtHashFwdIteratorGeneric<undefined_rank>,
                           QgramNtHashFwdIteratorGenericNoTransform
                             <undefined_rank>,
                           BinaryNtTable,
                           false>
                          (options.inputfile_get(),
                           options.qgram_length_get(),
                           options.s_get(),
                           options.r_get(),
                           options.num_threads_get());
    RunTimeClass rt_estimate{};
    const double F0 = table.estimate_F0();
    printf("F0\t%.0f\n", F0);
    printf("F1 (count)\t%zu\n",table.F1_count_get());
    printf("sequences_number\t%zu\n",table.sequences_number_get());
    rt_estimate.show("estimate");
  } else
  {
    NtTable table
      = is_protein ?
        ntcard_enumerate<split_at_wildcard,
                         QgramNtHashAAFwdIteratorGeneric<undefined_rank>,
                         QgramNtHashAAFwdIteratorGenericNoTransform
                           <undefined_rank>,
                         NtTable,
                         true>
                        (options.inputfile_get(),
                         options.qgram_length_get(),
                         options.s_get(),
                         options.r_get(),
                         options.num_threads_get())
        :
        ntcard_enumerate<split_at_wildcard,
                         QgramNtHashFwdIteratorGeneric<undefined_rank>,
                         QgramNtHashFwdIteratorGenericNoTransform
                           <undefined_rank>,
                         NtTable,
                         false>
                        (options.inputfile_get(),
                         options.qgram_length_get(),
                         options.s_get(),
                         options.r_get(),
                         options.num_threads_get());
    RunTimeClass rt_estimate{};
    if (options.fast_option_is_set())
    {
      printf("F0\t%.0f\n", table.estimate_F0());
      printf("F1 (count)\t%zu\n",table.F1_count_get());
    } else
    {
      const NtTableResult nt_table_result = table.estimate_all();
      if (options.show_f_option_is_set())
      {
        const size_t this_max = std::min(nt_table_result.t_max_get(),
                                         size_t(999));
        for (size_t idx = 1; idx <= this_max; idx++)
        {
          const double this_f_value = nt_table_result.f_n_at(idx);
          printf("f_%zu:\t%f\t%zu\n",
                  idx, this_f_value,
                  static_cast<size_t>(std::abs(this_f_value)));
        }
      }
      printf("F0\t%.0f\n",nt_table_result.F0_get());
      printf("F1 (estimate)\t%.0f\n",nt_table_result.F1_get());
      printf("F1 (count)\t%zu\n",nt_table_result.F1_count_get());
    }
    printf("sequences_number\t%zu\n",table.sequences_number_get());
    rt_estimate.show("estimate");
  }
}

int main(int argc, char **argv)
{
  RunTimeClass rt_all{};
  NtcardOptions options{};
  try
  {
    options.parse(argc, argv);
  }
  catch (const std::exception &err)
  {
    std::cerr << argv[0] << ": " << err.what() << std::endl;
    return EXIT_FAILURE;
  }
  if (options.help_option_is_set())
  {
    return EXIT_SUCCESS;
  }
  try
  {
    if (options.handle_wildcard_like_A_is_set())
    {
      estimate_F_values<false>(options);
    } else
    {
      estimate_F_values<true>(options);
    }
  }
  catch (const std::exception &err)
  {
    std::cerr << argv[0] << ": " << err.what() << std::endl;
    return EXIT_FAILURE;
  }
  StrFormat msg("ntcard.all\t%c\t%zu\t\t%s",
                options.binary_option_is_set() ? 'b' : 'n',
                options.num_threads_get(),
                options.inputfile_get().c_str());
  rt_all.show(msg.str());
  return EXIT_SUCCESS;
}

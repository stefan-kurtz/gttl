#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include "utilities/gttl_file_open.hpp"
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
  if (options.binary_option_is_set())
  {
    BinaryNtTable table
      = ntcard_enumerate<split_at_wildcard,
                         QgramNtHashFwdIteratorGeneric<undefined_rank>,
                         QgramNtHashFwdIteratorGenericNoTransform
                           <undefined_rank>,
                         BinaryNtTable>
                        (options.inputfile_get(),
                         options.qgram_length_get(),
                         options.s_get(),
                         options.r_get(),
                         options.num_threads_get());
    RunTimeClass rt_estimate{};
    const double F0 = table.estimate_F0();
    printf("F0 = %.0f\n", F0);
    rt_estimate.show("estimate");
  } else
  {
    NtTable table
      = ntcard_enumerate<split_at_wildcard,
                         QgramNtHashFwdIteratorGeneric<undefined_rank>,
                         QgramNtHashFwdIteratorGenericNoTransform
                           <undefined_rank>,
                         NtTable>
                        (options.inputfile_get(),
                         options.qgram_length_get(),
                         options.s_get(),
                         options.r_get(),
                         options.num_threads_get());
    RunTimeClass rt_estimate{};
    if (options.fast_option_is_set())
    {
      printf("F0 = %.0f\n", table.estimate_F0());
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
      printf("F0 = %.0f\n",nt_table_result.F0_get());
      printf("F1 (estimate) = %.0f\n",nt_table_result.F1_get());
      printf("F1 (count) = %zu\n",nt_table_result.F1_count_get());
    }
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
  catch (std::invalid_argument &e)
  {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
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
  catch (const std::string &msg)
  {
    std::cerr << argv[0] << ": " << msg << std::endl;
    return EXIT_FAILURE;
  }
  StrFormat msg("ntcard.all\t%c\t%zu\t\t%s",
                options.binary_option_is_set() ? 'b' : 'n',
                options.num_threads_get(),
                options.inputfile_get().c_str());
  rt_all.show(msg.str());
  return EXIT_SUCCESS;
}
#include <cassert>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cinttypes>
#include <cstdlib>
#include <exception>
#include <format>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include <utility>
#include "minimizer_opt.hpp"
#include "sequences/gttl_multiseq.hpp"
#include "sequences/qgrams_hash_nthash.hpp"
#include "sequences/gttl_minimizer_generator.hpp"
#include "utilities/runtime_class.hpp"

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
  throw std::runtime_error(
          std::format("cannot handle sequences_bits + hash_bits = {} + {} > 72",
                      sequences_bits,
                      requested_hash_bits));
}


template <class HashIterator, class MinimizerValueClass, class GeneratorT>
inline void run_minimizer_generator(const GttlMultiseq &multiseq,
                                    const MinimizerOptions &options,
                                    int hash_bits)
{
  const size_t qgram_length = options.qgram_length_get();
  const size_t window_size = options.window_size_get();
  const uint64_t hash_mask = (hash_bits == 60)
                               ? ~uint64_t{0}
                               : ((uint64_t{1} << hash_bits) - 1);

  const size_t nseq = multiseq.sequences_number_get();

  std::cout << "# Running GttlMinimizerGenerator test on " << nseq
            << " sequences\n";

  for (size_t seqid = 0; seqid < nseq; ++seqid)
  {
    const char* const seqptr = multiseq.sequence_ptr_get(seqid);
    const size_t seqlen = multiseq.sequence_length_get(seqid);

    std::cout << "# sequence " << seqid << " (len=" << seqlen << ")\n";

    GeneratorT gen(qgram_length, window_size, hash_mask, seqptr, seqlen, seqid);

    for (auto it : gen)
    {
      const MinimizerValueClass mv = *it;
      std::printf("%" PRIu64 "\t%zu\t%zu\n",
                  std::get<0>(mv),
                  std::get<1>(mv),
                  std::get<2>(mv));
    }
  }
}

int main(int argc, char *argv[])
{
  RunTimeClass rt{};
  MinimizerOptions options{};

  try
  {
    options.parse(argc, argv);
  }catch (const std::invalid_argument &err)
  {
    std::cerr << argv[0] << ": error " << err.what() << '\n';
    return EXIT_FAILURE;
  }

  if (options.help_option_is_set())
  {
    return EXIT_SUCCESS;
  }

  try
  {
    RunTimeClass rt_create_multiseq{};
    constexpr bool store_header = true;
    constexpr bool store_sequence = true;
    constexpr uint8_t padding_char = UINT8_MAX;

    const GttlMultiseq multiseq(options.inputfiles_get(),
                                store_header,
                                store_sequence,
                                padding_char,
                                options.canonical_option_is_set());
    rt_create_multiseq.show("reading input files and creating multiseq");

    for (auto &log : multiseq.statistics())
    {
      std::cout << "# " << log << '\n';
    }

    int hash_bits = -1;
    int var_sizeof_unit_hashed_qgram;
    assert(options.hash_bits_get() != -1);
    std::tie(hash_bits, var_sizeof_unit_hashed_qgram) =
        determine_hash_bits(multiseq.sequences_bits_get(),
                            options.hash_bits_get());

    assert(var_sizeof_unit_hashed_qgram == 8
           or var_sizeof_unit_hashed_qgram == 9);

    using MinimizerValueClass = std::tuple<uint64_t, size_t, size_t>;

    if (options.canonical_option_is_set())
    {
      using HashIter = QgramNtHashIterator4;
      using GeneratorT = GttlMinimizerGenerator<HashIter, MinimizerValueClass>;

      run_minimizer_generator<HashIter, MinimizerValueClass, GeneratorT>(
          multiseq, options, hash_bits);
    } else {
      using HashIter = QgramNtHashFwdIterator4;
      using GeneratorT = GttlMinimizerGenerator<HashIter, MinimizerValueClass>;

      run_minimizer_generator<HashIter, MinimizerValueClass, GeneratorT>(
          multiseq, options, hash_bits);
    }
  } catch (const std::exception &err)
  {
    for (auto &&inputfile : options.inputfiles_get())
    {
      std::cerr << argv[0] << ": file \"" << inputfile << "\" "
                << err.what() << '\n';
    }
    return EXIT_FAILURE;
  }

  rt.show("overall");
  return EXIT_SUCCESS;
}

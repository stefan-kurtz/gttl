#include <exception>
#include <cstdlib>
#include <cstddef>
#include <cassert>
#include <string>
#include <iostream>
#include <vector>
#include "sequences/char_finder.hpp"
#include "sequences/gttl_fasta_generator.hpp"
#include "sequences/guess_if_protein_seq.hpp"
#include "sequences/qgrams_hash_invint.hpp"
#include "sequences/unwords.hpp"
#include "unwords_opt.hpp"
#include "utilities/file_size.hpp"
#include "utilities/runtime_class.hpp"
#include "utilities/str_format.hpp"

int main(int argc, char *argv[])
{
  UnwordsOptions options;
  try
  {
    options.parse(argc, argv);
  }
  catch (const std::exception &err)
  {
    std::cerr << argv[0] << err.what() << '\n';
    return EXIT_FAILURE;
  }
  if (options.help_option_is_set())
  {
    return EXIT_SUCCESS;
  }
  bool guessed_protein_sequences = false;
  const std::vector<std::string> inputfiles = options.inputfiles_get();
  try
  {
    guessed_protein_sequences = guess_if_protein_file(inputfiles);
  }
  catch (const std::exception &err)
  {
    std::cerr << argv[0] << err.what() << '\n';
    return EXIT_FAILURE;
  }
  if (guessed_protein_sequences &&
      options.ignore_reverse_complement_option_is_set())
  {
    std::cerr << argv[0] << ": " << ("input seems to be aminoacid sequences, "
                                     "but option "
                                     "-i/--ignore_reverse_complement is only "
                                     "available for DNA sequences\n");
    return EXIT_FAILURE;
  }

  bool haserr = false;
  Unwords *unwords = nullptr;
  size_t qgram_length_max = 0;
  try
  {
    if (options.qgram_length_max_get() > 0)
    {
      qgram_length_max = options.qgram_length_max_get();
    } else
    {
      const size_t alphabetsize = guessed_protein_sequences ? size_t{20}
                                                            : size_t{4};
      size_t upperbound_sequence_length;
      if (options.ignore_reverse_complement_option_is_set())
      {
        assert(!guessed_protein_sequences);
        upperbound_sequence_length = gttl_file_size(inputfiles);
      } else
      {
        upperbound_sequence_length = 2 * gttl_file_size(inputfiles);
      }
      qgram_length_max = estimate_qgram_length_max(upperbound_sequence_length,
                                                   alphabetsize);
    }
    RunTimeClass rt_unwords_finder{};
    const int buf_size = size_t{1} << size_t{14};
    GttlFastAGenerator<buf_size> gttl_si(&inputfiles);
    if (options.store_sequences_option_is_set())
    {
      RunTimeClass rt_sequence_storing{};
      std::vector<GttlFastAEntry<buf_size>> sequences{};
      for (const auto *si : gttl_si)
      {
        sequences.emplace_back(std::string(""),
                               std::string(si->sequence_get()));
      }
      const StrFormat msg("storing %zu sequences", sequences.size());
      rt_sequence_storing.show(msg.str());
      unwords = unwords_finder<std::vector<GttlFastAEntry<buf_size>>>
                              (guessed_protein_sequences,
                               !options
                                 .ignore_reverse_complement_option_is_set(),
                               qgram_length_max,
                               sequences);
    } else
    {
      unwords = unwords_finder<GttlFastAGenerator<buf_size>>
                              (guessed_protein_sequences,
                               !options
                                 .ignore_reverse_complement_option_is_set(),
                               qgram_length_max,
                               gttl_si);
    }
    rt_unwords_finder.show("total");
  }
  catch (const std::exception &err)
  {
    std::cerr << argv[0] << ": " << err.what() << '\n';
    haserr = true;
  }
  if (!haserr)
  {
    if (unwords == nullptr)
    {
      std::cerr << argv[0] << ": all words of length <= " << qgram_length_max
                << " occur in the sequences, i.e. there are no unwords\n";
      haserr = true;
    } else
    {
      if (guessed_protein_sequences)
      {
        unwords->template show<20>(char_finder::aminoacids);
      } else
      {
        unwords->template
                 show<4>(InvertibleIntegercode2Iterator4::
                         alphabet.characters_get());
      }
    }
  }
  delete unwords;
  return haserr ? EXIT_FAILURE : EXIT_SUCCESS;
}

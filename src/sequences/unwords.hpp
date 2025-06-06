#ifndef UNWORDS_HPP
#define UNWORDS_HPP
#include <exception>
#include <vector>
#include <cmath>
#include "utilities/access_maybe_pointer.hpp"
#include "utilities/runtime_class.hpp"
#include "utilities/constexpr_for.hpp"
#include "utilities/str_format.hpp"
#include "utilities/multibitvector.hpp"
#include "sequences/char_range.hpp"
#include "sequences/char_finder.hpp"
#include "sequences/qgrams_hash_invint.hpp"
#include "sequences/qgram_decoder.hpp"

class Unwords
{
  size_t qgram_length,
         number_of_all_qgrams;
  Multibitvector<true> multibitvector;
  public:
  size_t sequences_total_length;
  Unwords(size_t alphabetsize,size_t _qgram_length)
    : qgram_length(_qgram_length)
    , number_of_all_qgrams(std::pow(alphabetsize,qgram_length))
    , multibitvector(number_of_all_qgrams)
    , sequences_total_length(0)
    {}
  [[nodiscard]] size_t qgram_length_get(void) const noexcept
  {
    return qgram_length;
  }
  [[nodiscard]] size_t number_of_all_qgrams_get(void) const noexcept
  {
    return number_of_all_qgrams;
  }
  [[nodiscard]] size_t size(void) const noexcept
  {
    return multibitvector.size() - multibitvector.count();
  }
  void set(uint64_t integer_code)
  {
    multibitvector.set(integer_code);
  }
  [[nodiscard]] std::vector<size_t> unwords_vector_get(void) const noexcept
  {
    std::vector<size_t> unwords_vector{};
    for (size_t integer_code = 0; integer_code < number_of_all_qgrams;
         integer_code++)
    {
      if (!multibitvector[integer_code])
      {
        unwords_vector.push_back(integer_code);
      }
    }
    assert(unwords_vector.size() == this->size());
    return unwords_vector;
  }
  template<size_t alphabetsize>
  void show(const char *alphabet) const noexcept
  {
    QgramDecoder<alphabetsize,false> qgram_decoder(alphabet,qgram_length);
    std::cout << "# alphabet size:\t" << alphabetsize << '\n';
    std::cout << "# length of qgrams:\t" << qgram_length << '\n';
    std::cout << "# number of all qgrams:\t" << number_of_all_qgrams
              << '\n';
    std::cout << "# total length of sequences:\t" << sequences_total_length
              << '\n';
    std::cout << "# number of unwords:\t" << size() << '\n';
    std::cout << "# unwords:\n";

    auto unwords_vector = unwords_vector_get();
    for (const auto integer_code : unwords_vector)
    {
      const char *qgram = qgram_decoder.decode(integer_code);
      std::cout << qgram << '\n';
    }
  }
};

template<class CharRanger, class InvertibleIntcodeIterator,
         bool reverse_complement_option, class SeqGenerator>
static Unwords *try_if_all_qgrams_occur(size_t qgram_length,
                                        size_t alphabetsize,
                                        SeqGenerator &seq_generator)
{
  Unwords *unwords = new Unwords(alphabetsize,qgram_length);
  bool all_qgrams_present = false;
  try
  {
    if constexpr(requires { seq_generator.reset(); })
    {
      seq_generator.reset();
    }
    for (auto const& si : seq_generator)
    {
      auto sequence = access_maybe_pointer(si).sequence_get();
      CharRanger ranger(sequence.data(),sequence.size());
      for (auto const &&range : ranger)
      {
        const size_t this_length = std::get<1>(range);
        unwords->sequences_total_length += this_length;
        const char *substring = sequence.data() + std::get<0>(range);
        InvertibleIntcodeIterator qgiter(qgram_length, substring,this_length);
        for (auto const &&code_pair : qgiter)
        {
          const uint64_t integer_code = std::get<0>(code_pair);
          unwords->set(integer_code);
          if constexpr (reverse_complement_option)
          {
            const uint64_t compl_integer_code = std::get<1>(code_pair);
            unwords->set(compl_integer_code);
          }
          if (unwords->size() == 0)
          {
            all_qgrams_present = true;
            break;
          }
        }
        if (all_qgrams_present)
        {
          break;
        }
      }
    }
  }
  catch (const std::exception &err)
  {
    delete unwords;
    throw;
  }
  return unwords;
}

template<class CharRanger, class InvertibleIntcodeIterator,
         bool reverse_complement, class SeqGenerator>
static Unwords *unwords_binary_search(size_t qgram_length_max,
                                      size_t alphabetsize,
                                      SeqGenerator &seq_generator)
{
  Unwords *last_successful_unwords = nullptr;
  size_t l = 1,
         r = qgram_length_max;
  while (l <= r)
  {
    const size_t qgram_length = (l+r)/2;

    RunTimeClass compute_unwords_runtime{};
    Unwords *unwords = try_if_all_qgrams_occur<CharRanger,
                                               InvertibleIntcodeIterator,
                                               reverse_complement,
                                               SeqGenerator>
                                              (qgram_length,
                                               alphabetsize,
                                               seq_generator);
    StrFormat msg("count number of different %zu-grams", qgram_length);
    compute_unwords_runtime.show(msg.str());
    if (unwords->size() > 0)
    {
      delete last_successful_unwords;
      last_successful_unwords = unwords;
      r = qgram_length-1;
    } else
    {
      l = qgram_length+1;
      delete unwords;
    }
  }
  return last_successful_unwords;
}

static constexpr const char_finder::NucleotideFinder unw_nucleotide_finder{};
static constexpr const char_finder::AminoacidFinder unw_aminoacid_finder{};

template<class SeqGenerator>
static Unwords *unwords_finder(bool is_protein_sequence,
                               bool reverse_complement,
                               size_t qgram_length_max,
                               SeqGenerator &seq_generator)
{
  Unwords *unwords = nullptr;
  if (!is_protein_sequence)
  {
    constexpr_for<0,1+1,1>([&](auto compile_time_reverse_complement)
    {
      if (compile_time_reverse_complement
          == static_cast<int>(reverse_complement))
      {
        using NucleotideRanger = GttlCharRange<char_finder::NucleotideFinder,
                                               unw_nucleotide_finder,
                                               true, false>;
        unwords = unwords_binary_search<NucleotideRanger,
                                        InvertibleIntegercode2Iterator4,
                                        compile_time_reverse_complement,
                                        SeqGenerator>
                                       (qgram_length_max,
                                        4,
                                        seq_generator);
      }
    });
  } else
  {
    using AminoacidRanger = GttlCharRange<char_finder::AminoacidFinder,
                                          unw_aminoacid_finder,
                                          true, false>;
    unwords = unwords_binary_search<AminoacidRanger,
                                    InvertibleIntegercodeIterator20,
                                    false,
                                    SeqGenerator>
                                   (qgram_length_max,
                                    20,
                                    seq_generator);
  }
  return unwords;
}

static size_t estimate_qgram_length_max(size_t upper_bound_sequence_length,
                                        size_t alphabetsize)
{
  const size_t qgram_length
    = std::ceil(std::log(static_cast<double>(upper_bound_sequence_length+1))/
                std::log(static_cast<double>(alphabetsize)));
  return qgram_length;
}

#ifdef EXACT_COMPUTATION_OF_SEQUENCE_LENGTH
template<bool reverse_complement, uint8_t alphabetsize>
size_t estimate_qgram_length_max(const std::vector<std::string> &inputfiles)
{
  size_t sequences_length = 0;
  const int buf_size = 1 << 14;
  GttlFastAGenerator<buf_size> gttl_si(&inputfiles);
  try
  {
    for (auto &&si : gttl_si)
    {
      sequences_length += si.sequence_get().length();
    }
  }
  catch (const std::exception &err)
  {
    throw;
  }
  if constexpr (reverse_complement)
  {
    sequences_length *= 2;
  }
  size_t qgram_length = std::ceil(
                        std::log(static_cast<double>(sequences_length+1))/
                        std::log(static_cast<double>(alphabetsize)));
  std::cout << "# total length of sequences:\t" << sequences_length
            << std::endl;
  std::cout << "# estimated length of unwords:\t" << qgram_length << std::endl;
  return qgram_length;
}
#endif
#endif

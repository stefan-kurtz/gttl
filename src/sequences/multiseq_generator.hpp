#ifndef MULTISEQ_GENERATOR_HPP
#define MULTISEQ_GENERATOR_HPP
#include <cstddef>
#include <cassert>
#include <utility>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>
#include "sequences/gttl_multiseq.hpp"

template <class GeneratorClass, bool Condition>
class OptionalGenerator;

template <class GeneratorClass>
class OptionalGenerator<GeneratorClass, true>
{
  GeneratorClass gen;
  GeneratorClass::Iterator it;
  public:
  OptionalGenerator(const char *inputfile)
    : gen(inputfile)
    , it(gen.begin())
  { }
  bool end(void)
  {
    return it == gen.end();
  }
  auto get(void) const
  {
    return *it;
  }
  void operator ++ (void)
  {
    ++it;
  }
};

template <class GeneratorClass>
class OptionalGenerator<GeneratorClass, false>
{
  public:
  OptionalGenerator([[maybe_unused]] const char *inputfile)
  { }
};

template<class SequenceGeneratorClass, bool fastq_paired_input>
class GttlMultiseqGenerator
{
  private:
  static constexpr const bool store_sequence = true;
  const std::vector<std::string> &inputfiles;
  GttlMultiseq *multiseq;
  OptionalGenerator<SequenceGeneratorClass, true> fastq_gen0;
  OptionalGenerator<SequenceGeneratorClass, fastq_paired_input> fastq_gen1;
  const size_t max_num_sequences;
  size_t sequence_number_offset;
  const uint8_t padding_char;
  const bool store_header;
  const bool has_owner_ship;
  bool plan_for_exhausted;
  const std::string unequal_length_error(const std::string &first_file,
                                         const std::string &second_file)
  {
    return std::string("paired readfiles have different number of reads: ") +
           first_file + std::string(" has less reads than ") + second_file;
  }
  public:
  explicit GttlMultiseqGenerator(const std::vector<std::string> &_inputfiles,
                                 bool _store_header,
                                 size_t _max_num_sequences,
                                 uint8_t _padding_char,
                                 bool grant_owner_ship)
    : inputfiles(_inputfiles)
    , multiseq(nullptr)
    , fastq_gen0(inputfiles[0].c_str())
    , fastq_gen1(fastq_paired_input ? inputfiles[1].c_str() : nullptr)
    , max_num_sequences(_max_num_sequences)
    , sequence_number_offset(0)
    , padding_char(_padding_char)
    , store_header(_store_header)
    , has_owner_ship(grant_owner_ship)
    , plan_for_exhausted(false)
  { }

  ~GttlMultiseqGenerator(void)
  {
    if (has_owner_ship and multiseq != nullptr)
    {
      delete multiseq;
      multiseq = nullptr;
    }
  }

  GttlMultiseq *data_get(void) const
  {
    return multiseq;
  }

  bool advance(void)
  {
    if (plan_for_exhausted)
    {
      return false;
    }
    if (has_owner_ship and multiseq != nullptr)
    {
      delete multiseq;
      multiseq = nullptr;
    }
    constexpr const bool has_reverse_complement = false;
    multiseq = new GttlMultiseq(store_sequence, /*CONSTRUCTOR */
                                padding_char,
                                sequence_number_offset,
                                fastq_paired_input,
                                has_reverse_complement);
    const size_t seqnum_upper_bound
      = fastq_paired_input ? (max_num_sequences/2)
                           : max_num_sequences;
    size_t seqnum;
    for (seqnum = 0; seqnum < seqnum_upper_bound; ++seqnum)
    {
      if (fastq_gen0.end())
      {
        if constexpr (fastq_paired_input)
        {
          if (not fastq_gen1.end())
          {
            throw std::runtime_error(unequal_length_error(inputfiles[0],
                                                          inputfiles[1]));
          }
        }
        plan_for_exhausted = true;
        break;
      }
      if constexpr (fastq_paired_input)
      {
        if (fastq_gen1.end())
        {
          throw std::runtime_error(unequal_length_error(inputfiles[1],
                                                        inputfiles[0]));
        }
      }
      multiseq->append(fastq_gen0.get()->header_get(),
                       fastq_gen0.get()->sequence_get(),
                       store_header,
                       store_sequence,
                       padding_char);
      ++fastq_gen0;
      if constexpr (fastq_paired_input)
      {
        multiseq->append(fastq_gen1.get()->header_get(),
                         fastq_gen1.get()->sequence_get(),
                         store_header,
                         store_sequence,
                         padding_char);
        ++fastq_gen1;
      }
    }
    if constexpr (fastq_paired_input)
    {
      sequence_number_offset += (2 * seqnum);
    } else
    {
      sequence_number_offset += seqnum;
    }
    return true;
  }
  class Iterator
  {
    private:
    GttlMultiseqGenerator *generator;
    bool iter_exhausted;
    public:
    explicit Iterator(GttlMultiseqGenerator* _generator, bool _exhausted)
      : generator(_generator)
      , iter_exhausted(_exhausted)
    {
      if (not iter_exhausted)
      {
        ++(*this);
      }
    }

    auto operator * (void) const
    {
      return generator->data_get();
    }

    const Iterator& operator ++ (void)
    {
      assert(not iter_exhausted);
      iter_exhausted = not generator->advance();
      return *this;
    }

    bool operator == (const Iterator& other) const
    {
      return iter_exhausted == other.iter_exhausted and
             generator == other.generator;
    }

    bool operator != (const Iterator& other) const
    {
      return not (*this == other);
    }
  };
  Iterator begin(void)
  {
    return Iterator(this, false);
  }

  Iterator end(void)
  {
    return Iterator(this, true);
  }
};
#endif

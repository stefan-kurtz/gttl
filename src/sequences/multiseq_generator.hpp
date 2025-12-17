#ifndef MULTISEQ_GENERATOR_HPP
#define MULTISEQ_GENERATOR_HPP
#include <cstddef>
#include <cassert>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>
#include "sequences/gttl_multiseq.hpp"
#include "utilities/string_values_join.hpp"

template <class GeneratorClass, bool Condition>
class OptionalGenerator;

template <class GeneratorClass>
class OptionalGenerator<GeneratorClass, true>
{
  GeneratorClass gen;
  GeneratorClass::Iterator it;
  public:
  OptionalGenerator(const std::vector<std::string> *inputfiles)
    : gen(inputfiles)
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
  OptionalGenerator([[maybe_unused]] const std::vector<std::string> *inputfiles)
  { }
};

template<class SequenceGeneratorClass, bool fastq_paired_input>
class GttlMultiseqGenerator
{
  private:
  using StrVec = std::vector<std::string>;
  [[nodiscard]]
  StrVec extract_file_list(size_t remainder, const StrVec &svec) const
  {
    assert(remainder == 0 or remainder == 1);
    StrVec selection;
    for (size_t idx = 0; idx < svec.size(); idx++)
    {
      if (idx % 2 == remainder)
      {
        selection.push_back(svec[idx]);
      }
    }
    return selection;
  }
  static constexpr const bool store_sequence = true;
  GttlMultiseq *multiseq;
  const StrVec file_list0;
  const StrVec file_list1;
  OptionalGenerator<SequenceGeneratorClass, true> gen0;
  OptionalGenerator<SequenceGeneratorClass, fastq_paired_input> gen1;
  const size_t max_num_sequences;
  size_t sequence_number_offset;
  const uint8_t padding_char;
  const bool store_header;
  const bool has_owner_ship;
  bool plan_for_exhausted;
  [[nodiscard]]
  const std::string unequal_length_error(const StrVec &f0, const StrVec &f1)
    const
  {
    const std::string fst_list = string_values_join(", ",f0.begin(),f0.end());
    const std::string snd_list = string_values_join(", ",f1.begin(),f1.end());
    return std::string("the read files ") + fst_list +
           std::string("have less reads than the read files ") + snd_list;
  }
  public:
  explicit GttlMultiseqGenerator(const StrVec &inputfiles,
                                 bool _store_header,
                                 size_t _max_num_sequences,
                                 uint8_t _padding_char,
                                 bool grant_owner_ship)
    : multiseq(nullptr)
    , file_list0(fastq_paired_input ? extract_file_list(size_t(0),inputfiles)
                                    : inputfiles)
    , file_list1(fastq_paired_input ? extract_file_list(size_t(1),inputfiles)
                                    : StrVec{})
    , gen0(&file_list0)
    , gen1(file_list1.empty() ? nullptr : &file_list1)
    , max_num_sequences(_max_num_sequences)
    , sequence_number_offset(0)
    , padding_char(_padding_char)
    , store_header(_store_header)
    , has_owner_ship(grant_owner_ship)
    , plan_for_exhausted(false)
  {
    if constexpr (fastq_paired_input)
    {
      if (file_list0.size() != file_list1.size())
      {
        throw std::invalid_argument("if the input consist of paired read files,"
                                    " the number of files must be even");
      }
    }
  }

  ~GttlMultiseqGenerator(void)
  {
    if (has_owner_ship and multiseq != nullptr)
    {
      delete multiseq;
      multiseq = nullptr;
    }
  }

  [[nodiscard]] GttlMultiseq *data_get(void) const
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
      if (gen0.end())
      {
        if constexpr (fastq_paired_input)
        {
          if (not gen1.end())
          {
            throw std::runtime_error(unequal_length_error(file_list0,
                                                          file_list1));
          }
        }
        plan_for_exhausted = true;
        break;
      }
      if constexpr (fastq_paired_input)
      {
        if (gen1.end())
        {
          throw std::runtime_error(unequal_length_error(file_list1,file_list0));
        }
      }
      multiseq->append(gen0.get()->header_get(),
                       gen0.get()->sequence_get(),
                       store_header,
                       store_sequence,
                       padding_char);
      ++gen0;
      if constexpr (fastq_paired_input)
      {
        multiseq->append(gen1.get()->header_get(),
                         gen1.get()->sequence_get(),
                         store_header,
                         store_sequence,
                         padding_char);
        ++gen1;
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
    GttlMultiseqGenerator *multiseq_generator;
    bool iter_exhausted;
    public:
    explicit Iterator(GttlMultiseqGenerator* _multiseq_generator,
                      bool _exhausted)
      : multiseq_generator(_multiseq_generator)
      , iter_exhausted(_exhausted)
    {
      if (not iter_exhausted)
      {
        ++(*this);
      }
    }

    auto operator * (void) const
    {
      return multiseq_generator->data_get();
    }

    const Iterator& operator ++ (void)
    {
      assert(not iter_exhausted);
      iter_exhausted = not multiseq_generator->advance();
      return *this;
    }

    bool operator == (const Iterator& other) const
    {
      return iter_exhausted == other.iter_exhausted and
             multiseq_generator == other.multiseq_generator;
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

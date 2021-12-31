#ifndef GTTL_MULTISEQ_HPP
#define GTTL_MULTISEQ_HPP

#include <tuple>
#include <iostream>
#include <string>
#include <string>
#include <string_view>
#include <vector>
#include <map>
#include <tuple>
#include <algorithm>
#include <climits>

#include "utilities/unused.hpp"
#include "utilities/mathsupport.hpp"
#include "sequences/gttl_seq_iterator.hpp"
#include "sequences/gttl_fastq_iterator.hpp"

/* A class to store various sequences and their header information.
 the inputfile is read in using fasta_reader.
 Constructor may throw std::string
 - std::range_error */

class GttlMultiseq
{
  private:

  std::vector<size_t> sequence_offsets{};
  std::string concatenated_sequences{};
  std::vector<std::string> headers{};
  size_t sequences_number{0},
         sequences_total_length{0},
         sequences_minimum_length{ULONG_MAX},
         sequences_maximum_length{0};
  uint8_t padding_char;
  bool constant_padding_char;
  std::map<size_t,size_t> length_dist_map{};

  public:

  template<bool store>
  void append(const std::string_view header,
              const std::string_view sequence,
              GTTL_UNUSED /* if !store */ uint8_t this_padding_char)
  {
    if constexpr (store)
    {
      headers.push_back(std::string(header.substr(1,header.size()-1)));
      concatenated_sequences += sequence;
      concatenated_sequences.push_back(static_cast<char>(this_padding_char));
      sequence_offsets.push_back(concatenated_sequences.size());
    }
    sequences_number++;
    sequences_total_length += sequence.size();
    sequences_minimum_length = std::min(sequence.size(),
                                        sequences_minimum_length);
    sequences_maximum_length = std::max(sequence.size(),
                                        sequences_maximum_length);
    length_dist_map[sequence.size()]++;
  }

  void append_padding_char(uint8_t this_padding_char)
  {
    concatenated_sequences.push_back(static_cast<char>(this_padding_char));
  }

  private:

  void init(uint8_t this_padding_char)
  {
    concatenated_sequences.push_back(static_cast<char>(this_padding_char));
    sequence_offsets.push_back(size_t(1));
  }

  void multiseq_reader(const std::vector<std::string> &inputfiles,
                       bool store, bool zip_readpair_files)
  {
    if (padding_char <= uint8_t(122))
    {
      throw std::string("non letter padding character exhausted");
    }
    if (store)
    {
      init(padding_char);
    }
    static constexpr const int buf_size = 1 << 14;
    if (zip_readpair_files)
    {
      assert(inputfiles.size() == 2);
      GttlFastQIterator<buf_size> fastq_it0(inputfiles[0]),
                                  fastq_it1(inputfiles[1]);

      auto it0 = fastq_it0.begin();
      auto it1 = fastq_it1.begin();
      if (store)
      {
        while (it0 != fastq_it0.end() && it1 != fastq_it1.end())
        {
          append<true>((*it0).header_get(),(*it0).sequence_get(),padding_char);
          append<true>((*it1).header_get(),(*it1).sequence_get(),padding_char);
          ++it0;
          ++it1;
        }
      } else
      {
        while (it0 != fastq_it0.end() && it1 != fastq_it1.end())
        {
          append<false>((*it0).header_get(),(*it0).sequence_get(),padding_char);
          append<false>((*it1).header_get(),(*it1).sequence_get(),padding_char);
          ++it0;
          ++it1;
        }
      }
      const bool fst_more = (it0 != fastq_it0.end() && it1 == fastq_it1.end());
      const bool snd_more = (it0 == fastq_it0.end() && it1 != fastq_it1.end());
      if (fst_more || snd_more)
      {
        StrFormat msg("processing readpair files %s and %s: %s file contains "
                      "more sequences than %s file",
                      inputfiles[0].c_str(),inputfiles[1].c_str(),
                      fst_more ? "first" : "second",
                      fst_more ? "second" : "fist");
        throw msg.str();
      }
    } else
    {
      GttlSeqIterator<buf_size> gttl_si(&inputfiles);
      if (store)
      {
        for (auto &&si : gttl_si)
        {
          append<true>(si.header_get(),si.sequence_get(),padding_char);
        }
      } else
      {
        for (auto &&si : gttl_si)
        {
          append<false>(si.header_get(),si.sequence_get(),padding_char);
        }
      }
    }
  }

  public:

  GttlMultiseq(const char *inputfile, bool store, uint8_t _padding_char)
    : padding_char(_padding_char),
      constant_padding_char(true)
  {
    std::vector<std::string> inputfiles{std::string(inputfile)};
    multiseq_reader(inputfiles,store,false);
  }

  GttlMultiseq(const std::string &inputfile, bool store, uint8_t _padding_char)
    : padding_char(_padding_char),
      constant_padding_char(true)
  {
    std::vector<std::string> inputfiles{inputfile};
    multiseq_reader(inputfiles,store,false);
  }

  GttlMultiseq(const std::vector<std::string> &inputfiles,bool store,
               uint8_t _padding_char)
    : padding_char(_padding_char),
      constant_padding_char(true)
  {
    multiseq_reader(inputfiles,store,false);
  }

  GttlMultiseq(const std::string &readpair_file1,
               const std::string &readpair_file2,
               bool store, uint8_t _padding_char)
    : padding_char(_padding_char),
      constant_padding_char(true)
  {
    std::vector<std::string> inputfiles{readpair_file1,readpair_file2};
    multiseq_reader(inputfiles,store,true);
  }

  GttlMultiseq(bool store, uint8_t _padding_char)
    : padding_char(_padding_char),
      constant_padding_char(true)
  {
    if (padding_char <= uint8_t(122))
    {
      throw std::string("non letter padding character exhausted");
    }
    if (store)
    {
      init(padding_char);
    }
  }

  GttlMultiseq(const std::vector<std::string> &inputfiles,
               const char *forbidden_as_padding)
    : constant_padding_char(false)
  {
    constexpr const uint8_t initial_padding_char = 0;
    init(initial_padding_char);
    uint32_t padding_char = initial_padding_char;
    bool mark[UINT8_MAX+1] = {false};
    for (const char *f = forbidden_as_padding; *f != '\0'; f++)
    {
      mark[static_cast<int>(*f)] = true;
    }
    static constexpr const int buf_size = 1 << 14;
    GttlSeqIterator<buf_size> gttl_si(&inputfiles);
    for (auto &&si : gttl_si)
    {
      assert(padding_char <= UINT8_MAX);
      append<true>(si.header_get(),si.sequence_get(),
                   static_cast<uint8_t>(padding_char));
      padding_char = (padding_char + 1) & uint32_t(255);
      padding_char += mark[padding_char];
    }
    /* in case the computation of the lcp is a suffix of the last sequence,
       so that the lcp-computation crosses sequence boundaries involving
       identical sequence padding characters. Then the position after
       the padding character for the last sequence is accessed. To provide
       a value, we add an additional padding char */
    append_padding_char(static_cast<uint8_t>(padding_char));
  }

  size_t sequences_number_get(void) const noexcept
  {
    return sequences_number;
  }

  size_t sequences_total_length_get(void) const noexcept
  {
    return sequences_total_length;
  }

  size_t sequences_minimum_length_get(void) const noexcept
  {
    return sequences_minimum_length;
  }

  size_t sequences_maximum_length_get(void) const noexcept
  {
    return sequences_maximum_length;
  }

  int sequences_length_bits_get(void) const noexcept
  {
    return gttl_required_bits<size_t>(sequences_maximum_length);
  }

  int sequences_number_bits_get(void) const noexcept
  {
    assert(sequences_number_get() > 0);
    return gttl_required_bits<size_t>(sequences_number_get() - 1);
  }

  int sequences_bits_get(void) const noexcept
  {
    return sequences_number_bits_get() + sequences_length_bits_get();
  }

  char padding_char_get(void) const
  {
    if (!constant_padding_char)
    {
      std::cerr << "programming error: " << __func__
                << " only works if the padding character is constant"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    return padding_char;
  }

  /* Give the length of sequence seqnum EXCLUDING padding symbol at the end */
  size_t sequence_length_get(size_t seqnum) const noexcept
  {
    /* To Check whether there are any problems considering the pointer at
    sequences_number goes out of bound and is used for the length of the last
    sequence */
    assert(seqnum + 1 < sequence_offsets.size() &&
           sequence_offsets[seqnum + 1] > sequence_offsets[seqnum]);
    return sequence_offsets[seqnum + 1] - sequence_offsets[seqnum] - 1;
  }

  /* Returns a pointer to the sequence with number seqnum */
  const char *sequence_ptr_get(size_t seqnum) const noexcept
  {
    assert(seqnum < sequences_number_get() &&
           sequence_offsets[seqnum] < concatenated_sequences.size());
    return concatenated_sequences.data() + sequence_offsets[seqnum];
  }

  const char *sequence_ptr_get(void) const noexcept
  {
    return sequence_ptr_get(0);
  }

  uint8_t sequence_char_get(size_t position) const noexcept
  {
    return static_cast<uint8_t>(concatenated_sequences[position+1]);
  }


  char *sequence_ptr_writable_get(size_t seqnum)
  {
    assert(seqnum < sequences_number_get() &&
           sequence_offsets[seqnum] < concatenated_sequences.size());
    return concatenated_sequences.data() + sequence_offsets[seqnum];
  }

  /* Returns length of header,
   * short version from start to first space(excluded). */
  size_t short_header_length_get(size_t seqnum) const noexcept
  {
    assert(seqnum < headers.size());
    size_t idx;
    for (idx = 0; idx < headers[seqnum].size() &&
                  !isspace(headers[seqnum][idx]); idx++)
      /* Nothing */;
    return idx;
  }

  const std::string_view header_get(size_t seqnum) const noexcept
  {
    assert(seqnum < headers.size());
    return headers[seqnum];
  }

  void statistics(void) const noexcept
  {
    std::cout << "# sequences_number\t" << sequences_number_get()
              << std::endl;
    std::cout << "# sequences_number_bits\t" << sequences_number_bits_get()
              << std::endl;
    std::cout << "# sequences_minimum_length\t"
              << sequences_minimum_length_get()
              << std::endl;
    std::cout << "# sequences_maximum_length\t"
              << sequences_maximum_length_get()
              << std::endl;
    std::cout << "# sequences_length_bits\t" << sequences_length_bits_get()
              << std::endl;
    std::cout << "# sequences_total_length\t" << sequences_total_length_get()
              << std::endl;
  }

  std::vector<std::pair<size_t,size_t>> length_distribution(void)
     const noexcept
  {
    std::vector<std::pair<size_t,size_t>> length_dist_table{};
    length_dist_table.reserve(length_dist_map.size());
    for (auto &&element : length_dist_map)
    {
      length_dist_table.push_back({std::get<0>(element),std::get<1>(element)});
    }
    std::sort(length_dist_table.begin(),length_dist_table.end());
    return length_dist_table;
  }
  size_t total_number_of_suffixes(size_t prefix_length) const noexcept
  {
    size_t total = 0;
    for (auto &&element : length_dist_map)
    {
      if (std::get<0>(element) >= prefix_length)
      {
        total += (std::get<0>(element) - prefix_length + 1) *
                 std::get<1>(element);
      }
    }
    return total;
  }
  /* Prints out the header and sequence infos to stdout
   - width gives the maximum line size, width 0
     prints prints out sequences in just one line
   - padding prints out the padding symbol if set to true,
   - small_header print out the header until the first space if set true
   - raw prints out sequences as saved in, raw=false converts ranks back to
     normal symbols for readability*/
  void show(size_t width, bool small_header) const noexcept
  {
    assert(concatenated_sequences.size() > 0);
    size_t seqnum;
#ifndef NDEBUG
    bool found_maximum_seq_length = false;
    bool found_minimum_seq_length = false;
#endif
    for (seqnum = 0; seqnum < sequences_number_get(); seqnum++)
    {
      std::cout.put('>');
      if (small_header)
      {
        size_t short_length = short_header_length_get(seqnum);
        std::cout << header_get(seqnum).substr(0,short_length);
      } else
      {
        std::cout << header_get(seqnum);
      }
      std::cout << std::endl;
      const char *currentseq = this->sequence_ptr_get(seqnum);
      size_t currentlength = this->sequence_length_get(seqnum);
      size_t line_width = 0;
      for (size_t idx = 0; idx < currentlength; idx++)
      {
        std::cout << currentseq[idx];
        line_width++;
        if (line_width == width)
        {
          std::cout << std::endl;
          line_width = 0;
        }
      }
      if (width == 0 || line_width > 0)
      {
        std::cout << std::endl;
      }
#ifndef NDEBUG
      assert(currentlength >= sequences_minimum_length);
      if (currentlength == sequences_minimum_length)
      {
        found_minimum_seq_length = true;
      }
      assert(currentlength <= sequences_maximum_length);
      if (currentlength == sequences_maximum_length)
      {
        found_maximum_seq_length = true;
      }
#endif
    }
    assert(found_minimum_seq_length);
    assert(found_maximum_seq_length);
  }

  /* Overload access operator[] */
  std::pair<const std::string_view,const std::string_view>
              operator[](size_t seqnum)
      const noexcept
  {
    assert(seqnum < sequences_number_get());
    const char *seq_ptr = sequence_ptr_get(seqnum);
    size_t seq_len = sequence_length_get(seqnum);
    std::string_view this_seq{seq_ptr,seq_len};

    return std::pair<std::string_view,std::string_view>(headers[seqnum],
                                                        this_seq);
  }
  template<class T,void (*transformation)(T &,char *,size_t)>
  void transformer(T &t)
  {
    for (size_t snum = 0; snum < sequences_number_get(); snum++)
    {
      transformation(t,sequence_ptr_writable_get(snum),
                       sequence_length_get(snum));
    }
  }
};

std::vector<GttlMultiseq *> multiseq_factory(const std::string &filename1,
                                             const std::string &filename2,
                                             size_t split_size)
{
  constexpr const int buf_size = 1 << 14;
  GttlFastQIterator<buf_size> fastq_it0(filename1),
                              fastq_it1(filename2);
  auto it0 = fastq_it0.begin();
  auto it1 = fastq_it1.begin();
  bool exhausted = false;
  std::vector<GttlMultiseq *> multiseq_vector{};

  while (!exhausted)
  {
    if (it0 == fastq_it0.end() || it1 == fastq_it1.end())
    {
      break;
    }
    const uint8_t padding_char = UINT8_MAX;
    GttlMultiseq *multiseq = new GttlMultiseq(true,padding_char);
    for (size_t idx = 0; idx < split_size; idx++)
    {
      if (it0 == fastq_it0.end() || it1 == fastq_it1.end())
      {
        exhausted = true;
        break;
      }
      multiseq->append<true>((*it0).header_get(),(*it0).sequence_get(),
                             padding_char);
      multiseq->append<true>((*it1).header_get(),(*it1).sequence_get(),
                             padding_char);
    }
    multiseq_vector.push_back(multiseq);
  }
  return multiseq_vector;
}
#endif  // MULTISEQ_HPP

#ifndef GTTL_MULTISEQ_HPP
#define GTTL_MULTISEQ_HPP

#include <cstdint>
#include <cassert>
#include <stdexcept>
#include <cctype>
#include <cstdlib>
#include <tuple>
#include <iostream>
#include <cstring>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <map>
#include <tuple>
#include <algorithm>
#include <climits>
#include <cstdio>

#include "sequences/gttl_fasta_generator.hpp"
#include "sequences/gttl_fastq_generator.hpp"
#include "utilities/str_format.hpp"
#include "utilities/mathsupport.hpp"
#include "utilities/cycle_of_numbers.hpp"

/* A class to store various sequences and their header information.
 the inputfile is read in using fasta_reader.
 Constructor may throw std::runtime_error
 - std::range_error */

class GttlMultiseq
{
  private:

  std::vector<size_t> sequence_offsets;
  std::string concatenated_sequences;
  std::vector<std::string> header_vector;
  size_t header_total_length{0},
         sequences_number{0},
         sequences_total_length{0},
         sequences_minimum_length{ULONG_MAX},
         sequences_maximum_length{0};
  uint8_t padding_char;
  bool constant_padding_char,
       has_reverse_complement,
       has_read_pairs;
  std::map<size_t,size_t> length_dist_map;
  std::vector<std::pair<uint16_t,uint16_t>> short_header_cache;

  void append_padding_char(uint8_t this_padding_char)
  {
    concatenated_sequences.push_back(static_cast<char>(this_padding_char));
  }

  public:

  void size_in_bytes_show(void) const noexcept
  {
    printf("Multiseq_size.class=%zu\n",sizeof(GttlMultiseq));
    printf("Multiseq_size.length_dist_map=%zu\n",
           length_dist_map.size() * 2 * sizeof(size_t));
    printf("Multiseq_size.header=%zu\n",
           header_vector.size() * sizeof(std::string) +
           header_total_length * sizeof(char));
    printf("Multiseq_size.seqoffset=%zu\n",
           sequence_offsets.size() * sizeof(size_t));
    printf("Multiseq_size.concatenated_sequences=%zu\n",
           concatenated_sequences.size() * sizeof(char));
  }

  [[nodiscard]] size_t size_in_bytes_extra(void) const noexcept
  {
    return sizeof(GttlMultiseq) +
           sizeof(std::string) * header_vector.size() +
           header_total_length * sizeof(char) +
           length_dist_map.size() * 2 * sizeof(size_t);
  }

  [[nodiscard]] size_t size_in_bytes_sequence(void) const noexcept
  {
    return sequence_offsets.size() * sizeof(size_t) +
           concatenated_sequences.size() * sizeof(char);
  }

  [[nodiscard]] size_t size_in_bytes(void) const noexcept
  {
    return size_in_bytes_sequence() + size_in_bytes_extra();
  }

  void append(const std::string_view header,
              const std::string_view sequence,
              bool store_header,
              bool store_sequence,
              [[maybe_unused]] /* if !store_sequence */
              uint8_t this_padding_char)
  {
    if (store_header)
    {
      header_vector.emplace_back(header);
      header_total_length += header.size();
    }
    if (store_sequence)
    {
      concatenated_sequences += sequence;
      append_padding_char(this_padding_char);
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

  private:

  void multipadding(bool add_offset,uint8_t this_padding_char, size_t how_many)
  {
    for (size_t idx = 0; idx < how_many; idx++)
    {
      append_padding_char(this_padding_char);
    }
    if (add_offset)
    {
      sequence_offsets.push_back(how_many);
    }
  }

  void padding_before_first_sequence(uint8_t this_padding_char)
  {
    multipadding(true,this_padding_char,7);
  }

  void multiseq_reader(const std::vector<std::string> &inputfiles,
                       bool store_header, bool store_sequence,
                       bool zip_readpair_files)
  {
    if (store_sequence)
    {
      padding_before_first_sequence(padding_char);
    }
    static constexpr const int buf_size = 1 << 14;
    if (zip_readpair_files)
    {
      assert(inputfiles.size() == 2);
      GttlFastQGenerator<buf_size> fastq_it0(inputfiles[0].c_str());
      GttlFastQGenerator<buf_size> fastq_it1(inputfiles[1].c_str());

      auto it0 = fastq_it0.begin();
      auto it1 = fastq_it1.begin();

      while (it0 != fastq_it0.end() && it1 != fastq_it1.end())
      {
        append((*it0)->header_get(), (*it0)->sequence_get(),
               store_header, store_sequence, padding_char);
        append((*it1)->header_get(),(*it1)->sequence_get(),
               store_header, store_sequence, padding_char);
        ++it0;
        ++it1;
      }
      const bool fst_more = (it0 != fastq_it0.end() && it1 == fastq_it1.end());
      const bool snd_more = (it0 == fastq_it0.end() && it1 != fastq_it1.end());
      if (fst_more || snd_more)
      {
        const StrFormat msg("processing readpair files %s and %s: %s file"
                            " contains more sequences than %s file",
                            inputfiles[0].c_str(),
                            inputfiles[1].c_str(),
                            fst_more ? "first" : "second",
                            fst_more ? "second" : "fist");
        throw std::runtime_error(msg.str());
      }
    } else
    {
      GttlFastAGenerator<buf_size> gttl_si(&inputfiles);
      for (auto &&si : gttl_si)
      {
        append(si->header_get(),si->sequence_get(),store_header,store_sequence,
               padding_char);
      }
    }
  }

  /* Returns start position and length of header substring,
   * short version from start to first space(excluded). */

  [[nodiscard]] std::pair<size_t, size_t>
  short_header_substring(const std::string_view header) const noexcept
  {
    size_t idx;
    for (idx = 0; idx < header.size() && !isspace(header[idx]); idx++)
      /* Nothing */ ;
    return std::make_pair(0,idx);
  }

  template <char first_delim, char second_delim>
  [[nodiscard]] std::pair<size_t, size_t>
  short_header_substring(const std::string_view header) const noexcept
  {
    const char *const first_delim_ptr = static_cast<const char *>(std::memchr(
                                 static_cast<const void *>(header.data()),
                                 first_delim,
                                 header.size()));
    if(first_delim_ptr != nullptr)
    {
      const char *const second_delim_ptr = static_cast<const char *>(memchr(
                                   static_cast<const void *>(first_delim_ptr
                                                             + 1),
                                   second_delim,
                                   static_cast<size_t>(header.data()
                                                       + header.size()
                                                       - (first_delim_ptr
                                                          + 1))));
      const char *const header_end       = second_delim_ptr != nullptr ?
                                         second_delim_ptr :
                                         (header.data() + header.size());
      return std::make_pair(static_cast<size_t>(first_delim_ptr + 1 -
                                                header.data()),
                            static_cast<size_t>(header_end -
                                                (first_delim_ptr+1)));
    }
    return short_header_substring(header);
  }

  public:

  void padding_after_last_sequence(uint8_t this_padding_char)
  {
    multipadding(false,this_padding_char,6);
  }

  GttlMultiseq(const char *inputfile, bool store_header, bool store_sequence,
               uint8_t _padding_char)
    : padding_char(_padding_char)
    , constant_padding_char(true)
    , has_reverse_complement(false)
    , has_read_pairs(false)
  {
    const std::vector<std::string> inputfiles{std::string(inputfile)};
    multiseq_reader(inputfiles,store_header, store_sequence, false);
  }

  GttlMultiseq(const std::string &inputfile,
               bool store_header, bool store_sequence,
               uint8_t _padding_char)
    : padding_char(_padding_char)
    , constant_padding_char(true)
    , has_reverse_complement(false)
    , has_read_pairs(false)
  {
    const std::vector<std::string> inputfiles{inputfile};
    multiseq_reader(inputfiles,store_header,store_sequence,false);
  }

  GttlMultiseq(const std::vector<std::string> &inputfiles,
               bool store_header, bool store_sequence,
               uint8_t _padding_char)
    : padding_char(_padding_char)
    , constant_padding_char(true)
    , has_reverse_complement(false)
    , has_read_pairs(false)
  {
    multiseq_reader(inputfiles,store_header,store_sequence,false);
  }

  GttlMultiseq(const std::string &readpair_file1,
               const std::string &readpair_file2,
               bool store_header, bool store_sequence, uint8_t _padding_char)
    : padding_char(_padding_char)
    , constant_padding_char(true)
    , has_reverse_complement(false)
    , has_read_pairs(true)
  {
    const std::vector<std::string> inputfiles{readpair_file1, readpair_file2};
    constexpr const bool zip_readpair_files = true;
    multiseq_reader(inputfiles,store_header,store_sequence,zip_readpair_files);
  }

  GttlMultiseq(bool store_sequence, uint8_t _padding_char)
    : padding_char(_padding_char)
    , constant_padding_char(true)
    , has_reverse_complement(false)
    , has_read_pairs(false)
  {
    if (store_sequence)
    {
      padding_before_first_sequence(padding_char);
    }
  }

  GttlMultiseq(const std::vector<std::string> &inputfiles,
               const std::vector<uint8_t> &forbidden_as_padding,
               bool store_header)
    : constant_padding_char(false)
    , has_reverse_complement(false)
    , has_read_pairs(false)
  {
    CycleOfNumbers cycle_of_numbers(forbidden_as_padding);
    uint8_t this_padding_char = cycle_of_numbers.next();
    padding_before_first_sequence(this_padding_char);
    static constexpr const int buf_size = 1 << 14;
    GttlFastAGenerator<buf_size> gttl_si(&inputfiles);
    for (auto &&si : gttl_si)
    {
      this_padding_char = cycle_of_numbers.next();
      constexpr const bool store_sequence = true;
      append(si->header_get(),si->sequence_get(),store_header,store_sequence,
             this_padding_char);
    }
    /* in case the computation of the lcp is a suffix of the last sequence,
       so that the lcp-computation crosses sequence boundaries involving
       identical sequence padding characters. Then the positions after
       the padding character for the last sequence are accessed. When
       blockwise comparisons are used, even 7 characters after the end
       of the last sequence are accessed. So after the padding character
       used anyway, we add 6 more padding characters. */
    padding_after_last_sequence(this_padding_char);
  }

  ~GttlMultiseq(void)
  { }

  [[nodiscard]] size_t sequences_number_get(void) const noexcept
  {
    return sequences_number;
  }

  [[nodiscard]] size_t sequences_total_length_get(void) const noexcept
  {
    return sequences_total_length;
  }

  [[nodiscard]] size_t sequences_minimum_length_get(void) const noexcept
  {
    return sequences_minimum_length;
  }

  [[nodiscard]] size_t sequences_maximum_length_get(void) const noexcept
  {
    return sequences_maximum_length;
  }

  [[nodiscard]] int sequences_length_bits_get(void) const noexcept
  {
    return gttl_required_bits<size_t>(sequences_maximum_length);
  }

  [[nodiscard]] int sequences_number_bits_get(void) const noexcept
  {
    assert(sequences_number_get() > 0);
    return gttl_required_bits<size_t>(sequences_number_get() - 1);
  }

  [[nodiscard]] int sequences_bits_get(void) const noexcept
  {
    return sequences_number_bits_get() + sequences_length_bits_get();
  }

  [[nodiscard]] char padding_char_get(void) const
  {
    if (!constant_padding_char)
    {
      std::cerr << "programming error: " << __func__
                << " only works if the padding character is constant\n";
      exit(EXIT_FAILURE);
    }
    return padding_char;
  }

  void has_reverse_complement_set(void) noexcept
  {
    has_reverse_complement = true;
  }

  [[nodiscard]] bool has_reverse_complement_is_set(void) const noexcept
  {
    return has_reverse_complement;
  }

  [[nodiscard]] bool has_read_pairs_is_set(void) const noexcept
  {
    return has_read_pairs;
  }

  /* Give the length of sequence seqnum EXCLUDING padding symbol at the end */
  [[nodiscard]] size_t sequence_length_get(size_t seqnum) const noexcept
  {
    /* To Check whether there are any problems considering the pointer at
    sequences_number goes out of bound and is used for the length of the last
    sequence */
    assert(seqnum + 1 < sequence_offsets.size() &&
            sequence_offsets[seqnum + 1] > sequence_offsets[seqnum]);
    return sequence_offsets[seqnum + 1] - sequence_offsets[seqnum] - 1;
  }

  /* Returns a pointer to the sequence with number seqnum */
  [[nodiscard]] const char *sequence_ptr_get(size_t seqnum) const noexcept
  {
    assert(seqnum < sequences_number_get() &&
            sequence_offsets[seqnum] < concatenated_sequences.size());
    return concatenated_sequences.data() + sequence_offsets[seqnum];
  }

  /* Returns a pointer to the sequence with number seqnum */
  /* This function shall only be called after transforming the sequences
      using a LiterateMultiseq */
  [[nodiscard]] const uint8_t *
  encoded_sequence_ptr_get(size_t seqnum) const noexcept
  {
    assert(seqnum < sequences_number_get() &&
            sequence_offsets[seqnum] < concatenated_sequences.size());
    return reinterpret_cast<const uint8_t *>
                            (concatenated_sequences.data() +
                            sequence_offsets[seqnum]);
  }

  [[nodiscard]] const char *sequence_ptr_get(void) const noexcept
  {
    return sequence_ptr_get(0);
  }

  [[nodiscard]] uint8_t sequence_char_get(size_t position) const noexcept
  {
    return static_cast<uint8_t>(concatenated_sequences[position+1]);
  }

  char *sequence_ptr_writable_get(size_t seqnum)
  {
    assert(seqnum < sequences_number_get() &&
            sequence_offsets[seqnum] < concatenated_sequences.size());
    return concatenated_sequences.data() + sequence_offsets[seqnum];
  }

  [[nodiscard]] const std::string_view header_get(size_t seqnum) const noexcept
  {
    assert(seqnum < header_vector.size());
    return header_vector[seqnum];
  }

  [[nodiscard]] std::pair<size_t, size_t>
  short_header_get(size_t seqnum) const noexcept
  {
    assert(seqnum < short_header_cache.size());
    uint16_t sh_offset;
    uint16_t sh_len;
    std::tie(sh_offset,sh_len) = short_header_cache[seqnum];
    return std::make_pair(static_cast<size_t>(sh_offset),
                          static_cast<size_t>(sh_len));
  }

  [[nodiscard]] std::vector<std::string> statistics() const noexcept
  {
    std::vector<std::string> log_vector;
    log_vector.push_back(std::string("sequences_number\t") +
                          std::to_string(sequences_number_get()));
    log_vector.push_back(std::string("sequences_number_bits\t") +
                          std::to_string(sequences_number_bits_get()));
    log_vector.push_back(std::string("sequences_minimum_length\t") +
                          std::to_string(sequences_minimum_length_get()));
    log_vector.push_back(std::string("sequences_maximum_length\t") +
                          std::to_string(sequences_maximum_length_get()));
    log_vector.push_back(std::string("sequences_length_bits\t") +
                          std::to_string(sequences_length_bits_get()));
    log_vector.push_back(std::string("sequences_total_length\t") +
                          std::to_string(sequences_total_length_get()));
    return log_vector;
  }

  [[nodiscard]] std::vector<std::pair<size_t, size_t>>
  length_distribution(size_t length_dist_bin_size = 1) const noexcept
  {
    std::vector<std::pair<size_t,size_t>> length_dist_table;
    if (length_dist_bin_size == 1)
    {
      length_dist_table.reserve(length_dist_map.size());
      for (auto const& [length,  count] : length_dist_map)
      {
        length_dist_table.push_back(std::make_pair(length,count));
      }
    } else
    {
      std::map<size_t,size_t> length_dist_bins;
      for (auto const& [length,  count] : length_dist_map)
      {
        length_dist_bins[length/length_dist_bin_size] += count;
      }
      length_dist_table.reserve(length_dist_bins.size());
      for (auto const& [length,  count] : length_dist_bins)
      {
        length_dist_table.push_back(std::make_pair(length,count));
      }
    }
    std::sort(length_dist_table.begin(),length_dist_table.end());
    return length_dist_table;
  }
  [[nodiscard]] size_t
  total_number_of_suffixes(size_t prefix_length) const noexcept
  {
    size_t total = 0;
    for (auto const& [length, count] : length_dist_map)
    {
      if (length >= prefix_length)
      {
        total += (length - prefix_length + 1) * count;
      }
    }
    return total;
  }
  [[nodiscard]] std::pair<const char *, size_t>
  header_ptr_with_length(size_t seqnum, bool short_header) const
  {
    const std::string_view header = header_get(seqnum);
    size_t header_offset;
    size_t header_len;
    if (short_header)
    {
      std::tie(header_offset,header_len) = short_header_get(seqnum);
    } else
    {
      header_offset = 0;
      header_len = header.size();
    }
    return std::make_pair(header.data() + header_offset,header_len);
  }
  void show_single_sequence(size_t width, bool short_header, size_t seqnum)
        const noexcept
  {
    printf(">");
    const char *header_ptr;
    size_t header_len;
    std::tie(header_ptr,header_len)
      = header_ptr_with_length(seqnum,short_header);
    std::fwrite(header_ptr,sizeof (char),header_len,stdout);
    printf("\n");
    const char *const currentseq = this->sequence_ptr_get(seqnum);
    const size_t currentlength = this->sequence_length_get(seqnum);
    size_t line_width = 0;
    for (size_t idx = 0; idx < currentlength; idx++)
    {
      printf("%c",currentseq[idx]);
      line_width++;
      if (line_width == width)
      {
        printf("\n");
        line_width = 0;
      }
    }
    if (width == 0 || line_width > 0)
    {
      printf("\n");
    }
  }
  /* Prints out the header and sequence infos to stdout
   - width gives the maximum line size, width 0
     prints prints out sequences in just one line
   - padding prints out the padding symbol if set to true,
   - short_header print out the header until the first space if set true
   - raw prints out sequences as saved in, raw=false converts ranks back to
     normal symbols for readability*/
  void show(size_t width, bool short_header) const noexcept
  {
    assert(concatenated_sequences.size() > 0);
#ifndef NDEBUG
    bool found_maximum_seq_length = false;
    bool found_minimum_seq_length = false;
#endif
    for (size_t seqnum = 0; seqnum < sequences_number_get(); seqnum++)
    {
      show_single_sequence(width,short_header,seqnum);
#ifndef NDEBUG
      const size_t currentlength = this->sequence_length_get(seqnum);
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
  void show_sorted_by_header(size_t width, bool short_header) const
  {
    assert(concatenated_sequences.size() > 0);
    std::vector<std::pair<std::string,size_t>> header_with_seqnum;
    header_with_seqnum.reserve(sequences_number_get());
    for (size_t seqnum = 0; seqnum < sequences_number_get(); seqnum++)
    {
      const char *header_ptr;
      size_t header_len;
      std::tie(header_ptr,header_len)
        = header_ptr_with_length(seqnum,short_header);
      const std::string header_substring(header_ptr, header_len);
      header_with_seqnum.push_back(std::make_pair(header_substring,seqnum));
    }
    std::sort(header_with_seqnum.begin(),header_with_seqnum.end());
    for (size_t idx = 1; idx < header_with_seqnum.size(); idx++)
    {
      assert(std::get<0>(header_with_seqnum[idx-1]) <=
             std::get<0>(header_with_seqnum[idx]));
      if (std::get<0>(header_with_seqnum[idx-1]) ==
          std::get<0>(header_with_seqnum[idx]))
      {
        throw std::runtime_error(std::string("sequence set contains a "
                                             "duplicated header ") +
                                 std::get<0>(header_with_seqnum[idx]));
      }
    }
    for (auto &&hws : header_with_seqnum)
    {
      show_single_sequence(width,short_header,std::get<1>(hws));
    }
  }

  /* Overload access operator[] */
  std::string_view operator [](size_t idx) const noexcept
  {
    const char *const seq_ptr = sequence_ptr_get(idx);
    const size_t len = sequence_length_get(idx);
    return std::string_view(seq_ptr,len);
  }

  template<class T,void (*transformation)(T &,char *,size_t)>
  void transformer(T &t)
  {
    for (size_t seqnum = 0; seqnum < sequences_number_get(); seqnum++)
    {
      transformation(t,sequence_ptr_writable_get(seqnum),
                       sequence_length_get(seqnum));
    }
  }

  template<char first_delim,char second_delim>
  void short_header_cache_create(void)
  {
    assert(sequences_number_get() > 0);
    for (size_t seqnum = 0; seqnum < sequences_number_get(); seqnum++)
    {
      const std::string_view header = header_get(seqnum);
      size_t sh_offset;
      size_t sh_len;
      std::tie(sh_offset,sh_len)
        = short_header_substring<first_delim,second_delim>(header);
      assert(sh_offset <= UINT16_MAX && sh_len <= UINT16_MAX);
      short_header_cache.push_back(std::make_pair(
                                     static_cast<uint16_t>(sh_offset),
                                     static_cast<uint16_t>(sh_len)));
    }
  }

  void short_header_cache_create(void)
  {
    assert(sequences_number_get() > 0);
    for (size_t seqnum = 0; seqnum < sequences_number_get(); seqnum++)
    {
      const std::string_view header = header_get(seqnum);
      size_t sh_offset;
      size_t sh_len;
      std::tie(sh_offset,sh_len)
        = short_header_substring(header);
      assert(sh_offset <= UINT16_MAX && sh_len <= UINT16_MAX);
      short_header_cache.push_back(std::make_pair(
                                     static_cast<uint16_t>(sh_offset),
                                     static_cast<uint16_t>(sh_len)));
    }
  }
};

template<char (*complement_base)(char)>
static inline GttlMultiseq *multiseq_with_reverse_complement(
                           const std::vector<std::string> &inputfiles,
                           bool store_header,
                           bool store_sequence,
                           uint8_t padding_char)
{
  static constexpr const int buf_size = 1 << 14;
  GttlFastAGenerator<buf_size> gttl_si(&inputfiles);
  GttlMultiseq *multiseq
    = new GttlMultiseq(store_sequence,padding_char);/*CONSTRUCTOR*/
  multiseq->has_reverse_complement_set();
  for (auto &&si : gttl_si)
  {
    const std::string_view &seq = si->sequence_get();
    multiseq->append(si->header_get(),seq,store_header,store_sequence,
                     padding_char);
    std::string rc_seq;
    rc_seq.reserve(seq.size());
    size_t bck = seq.size();
    while (bck > 0)
    {
      bck--;
      rc_seq.push_back(complement_base(seq[bck]));
    }
    multiseq->append(si->header_get(),rc_seq,store_header,store_sequence,
                     padding_char);
  }
  return multiseq;
}
#endif  // MULTISEQ_HPP

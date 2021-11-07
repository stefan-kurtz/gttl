#ifndef GTTL_MULTISEQ_HPP
#define GTTL_MULTISEQ_HPP

#include <tuple>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <climits>

#include "utilities/mathsupport.hpp"
#include "utilities/gttl_file_open.hpp"
#include "sequences/gttl_seq_iterator.hpp"

/* A class to store various sequences and their header information.
 the inputfile is read in using fasta_reader.
 Constructor may throw std::string
 - std::range_error */

class GttlMultiseq
{
  struct Iterator
  {
    private:
      char **sequence_ptr;
      size_t seqnum;
    public:
      Iterator(char **_sequence_ptr,size_t _seqnum) :
        sequence_ptr(_sequence_ptr),
        seqnum(_seqnum) {}
      std::string operator*(void) const
      {
        return std::string(sequence_ptr[seqnum],
                           static_cast<size_t>(sequence_ptr[seqnum + 1] -
                                               sequence_ptr[seqnum] - 1));
      }
      Iterator& operator++() /* prefix increment*/
      {
        seqnum++;
        return *this;
      }
      bool operator != (const Iterator& other) const noexcept
      {
        return seqnum != other.seqnum;
      }
  };
  private:
  bool store;
  size_t sequences_number,
         sequences_total_length,
         sequences_maximum_length,
         headers_total_length;
  char **sequence_ptr, **header_ptr;
  int sequences_number_bits,
      sequences_length_bits;

  uint8_t padding_char;

  void multiseq_reader(const std::vector<std::string> &inputfiles)
  {
    if (padding_char <= uint8_t(122))
    {
      throw std::string("non letter padding character exhausted");
    }
    static constexpr const int buf_size = (((size_t) 1) << 14);
    if (store)
    {
      for (auto && inputfile : inputfiles)
      {
        GttlFpType in_fp = gttl_fp_type_open(inputfile.c_str(), "rb");
        if (in_fp == nullptr)
        {
          throw std::string(": cannot open file");
        }
        GttlSeqIterator<buf_size> gttl_si_first_pass(in_fp);
        try
        {
          for (auto &&si : gttl_si_first_pass)
          {
            std::string header = std::get<0>(si);
            std::string sequence = std::get<1>(si);

            sequences_total_length += sequence.size();
            sequences_maximum_length = std::max(sequences_maximum_length,
                                                sequence.size());
            sequences_number++;
            assert(header.size() >= 2 && header[0] == '>' &&
                   header[header.size() - 1] == '\n');
            headers_total_length += (header.size() - 2);
          }
        }
        catch (std::string &msg)
        {
          gttl_fp_type_close(in_fp);
          throw msg;
        }
        gttl_fp_type_close(in_fp);
      }
      sequence_ptr
        = static_cast<char **>(malloc((sequences_number + 1) *
                                      sizeof *sequence_ptr));
      header_ptr
        = static_cast<char **>(malloc((sequences_number + 1) *
                                      sizeof *header_ptr));

      /* sequence_ptr[0] is malloced to area of size for all symbols
         (sequences_total_length) and padding sumbols, one per sequence,
         therefore sequences_number, plus one for the padding symbol at
         the start*/
      const size_t sizeof_mem_for_seqs
        = (sequences_total_length + sequences_number + 1) *
          sizeof *sequence_ptr;
      const int return_value
        = posix_memalign(reinterpret_cast<void **>(&sequence_ptr[0]), 32,
                                                   sizeof_mem_for_seqs);
      if (return_value == ENOMEM)
      {
        free(header_ptr);
        free(sequence_ptr);
        StrFormat msg("cannot allocate %lu bytes for all sequences",
                      sizeof_mem_for_seqs);
        throw msg.str();
      }
      sequence_ptr[0][0] = static_cast<char>(padding_char);
      sequence_ptr[0]++;
      header_ptr[0] = static_cast<char *>(malloc(headers_total_length *
                                                 sizeof *header_ptr));
      if (header_ptr[0] == NULL)
      {
        free(sequence_ptr[0]);
        free(sequence_ptr);
        free(header_ptr);
        StrFormat msg("cannot allocated %lu bytes for all headers",
                      headers_total_length);
        throw msg.str();
      }
      size_t seqnum = 0;
      for (auto && inputfile : inputfiles)
      {
        GttlFpType in_fp = gttl_fp_type_open(inputfile.c_str(), "rb");
        if (in_fp == nullptr)
        {
          throw std::string(": cannot open file");
        }
        GttlSeqIterator<buf_size> gttl_si_second_pass(in_fp);
        for (auto &&si : gttl_si_second_pass)
        {
          std::string header = std::get<0>(si);
          std::string sequence = std::get<1>(si);

          memcpy(sequence_ptr[seqnum], sequence.data(), sequence.size());
          *(sequence_ptr[seqnum] + sequence.size())
            = static_cast<char>(padding_char);
          sequence_ptr[seqnum + 1] = sequence_ptr[seqnum] + sequence.size() + 1;
          memcpy(header_ptr[seqnum], header.data() + 1, header.size() - 2);
          header_ptr[seqnum + 1] = header_ptr[seqnum] + header.size() - 2;
          seqnum++;
        }
        gttl_fp_type_close(in_fp);
      }
    } else
    {
      for (auto && inputfile : inputfiles)
      {
        GttlFpType in_fp = gttl_fp_type_open(inputfile.c_str(), "rb");
        if (in_fp == nullptr)
        {
          throw std::string(": cannot open file");
        }
        GttlSeqIterator<buf_size> gttl_si_first_pass(in_fp);
        try /* need this, as the catch needs to close the file pointer
               to prevent a memory leak */
        {
          for (auto &&si : gttl_si_first_pass)
          {
            auto sequence = std::get<1>(si);
            sequences_total_length += sequence.size();
            sequences_maximum_length = std::max(sequences_maximum_length,
                                                sequence.size());
            sequences_number++;
          }
        }
        catch (std::string &msg)
        {
          gttl_fp_type_close(in_fp);
          throw msg;
        }
        gttl_fp_type_close(in_fp);
      }
    }
    sequences_number_bits = gt_required_bits(sequences_number - 1);
    sequences_length_bits = gt_required_bits(sequences_maximum_length);
  }
  public:
  /* Constructor
   Inputfile should be in Fasta format, throws std::string,
   if multiple GttlMultseq-instance are used and there are pairwise
   comparisons of the sequences, use a different padding_char, so
   that one does not have to have a special case for handling sequence
   boundaries when for example computing maximal matches. */
  GttlMultiseq(const char *inputfile, bool _store, uint8_t _padding_char)
      : store(_store),
        sequences_number(0),
        sequences_total_length(0),
        sequences_maximum_length(0),
        headers_total_length(0),
        padding_char(_padding_char)
  {
    std::vector<std::string> inputfiles{};
    inputfiles.push_back(std::string(inputfile));
    multiseq_reader(inputfiles);
  }
  GttlMultiseq(const std::vector<std::string> &inputfiles,bool _store,
               uint8_t _padding_char)
      : store(_store),
        sequences_number(0),
        sequences_total_length(0),
        sequences_maximum_length(0),
        headers_total_length(0),
        padding_char(_padding_char)
  {
    multiseq_reader(inputfiles);
  }
  ~GttlMultiseq(void)
  {
    if (store)
    {
      sequence_ptr[0]--;
      free(sequence_ptr[0]);
      free(header_ptr[0]);
      free(sequence_ptr);
      free(header_ptr);
    }
  }

  size_t sequences_number_get(void) const noexcept
  {
    return sequences_number;
  }

  size_t sequences_total_length_get(void) const noexcept
  {
    return sequences_total_length;
  }

  size_t sequences_maximum_length_get(void) const
  {
    return sequences_maximum_length;
  }

  int sequences_length_bits_get(void) const noexcept
  {
    return sequences_length_bits;
  }

  int sequences_number_bits_get(void) const noexcept
  {
    return sequences_number_bits;
  }

  int sequences_bits_get(void) const noexcept
  {
    return sequences_number_bits_get() + sequences_length_bits_get();
  }

  uint8_t padding_char_get(void) const noexcept
  {
    assert(store);
    return padding_char;
  }

  /* Give the length of sequence seqnum EXCLUDING padding symbol at the end */
  size_t sequence_length_get(size_t seqnum) const noexcept
  {
    assert(store && seqnum < sequences_number_get() &&
           sequence_ptr[seqnum + 1] > sequence_ptr[seqnum]);
    /* To Check whether there are any problems considering the pointer at
    sequences_number goes out of bound and is used for the length of the last
    sequence */
    return static_cast<size_t>(sequence_ptr[seqnum + 1] -
                               sequence_ptr[seqnum] - 1);
  }

  /* Returns a pointer to the sequence with number seqnum */
  char *sequence_ptr_get(size_t seqnum) const noexcept
  {
    assert(store && seqnum < sequences_number_get());
    return sequence_ptr[seqnum];
  }

  char *sequence_ptr_get() const noexcept
  {
    return sequence_ptr[0];
  }

  /* Returns length of header. */
  size_t header_length_get(size_t seqnum) const noexcept
  {
    assert(store && seqnum < sequences_number_get() &&
           header_ptr[seqnum] < header_ptr[seqnum + 1]);
    return static_cast<size_t>(header_ptr[seqnum + 1] - header_ptr[seqnum]);
  }

  /* Returns length of header,
   * short version from start to first space(excluded). */
  size_t short_header_length_get(size_t seqnum) const noexcept
  {
    assert(store && seqnum < sequences_number_get() &&
           header_ptr[seqnum] < header_ptr[seqnum + 1]);
    const char *itr;
    for (itr = header_ptr[seqnum];
         itr < header_ptr[seqnum + 1] && !isspace(*itr); itr++)
      /* Nothing */;
    return static_cast<size_t>(itr - header_ptr[seqnum]);
  }

  /* Returns a pointer to the header of sequence number seqnum */
  const char *header_ptr_get(size_t seqnum) const noexcept
  {
    assert(store && seqnum < sequences_number_get());
    return header_ptr[seqnum];
  }

  void statistics(void)
  {
    std::cout << "# sequences_number\t" << sequences_number_get()
              << std::endl;
    std::cout << "# sequences_number_bits\t" << sequences_number_bits_get()
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
    std::map<size_t,size_t> length_dist_map{};

    for (size_t seqnum = 0; seqnum < sequences_number_get(); seqnum++)
    {
      size_t this_length = sequence_length_get(seqnum);
      if (length_dist_map.count(this_length) == 0)
      {
        length_dist_map[this_length] = 1;
      } else
      {
        length_dist_map[this_length]++;
      }
    }
    std::vector<std::pair<size_t,size_t>> length_dist_table{};
    length_dist_table.reserve(length_dist_map.size());
    for (auto &&element : length_dist_map)
    {
      length_dist_table.push_back({std::get<0>(element),std::get<1>(element)});
    }
    std::sort(length_dist_table.begin(),length_dist_table.end());
    return length_dist_table;
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
    assert(store);
    size_t seqnum;
#ifndef NDEBUG
    bool found_maximum_seq_length = false;
#endif
    for (seqnum = 0; seqnum < sequences_number_get(); seqnum++)
    {
      std::cout.put('>');
      if (small_header)
      {
        std::cout.write((const char *) this->header_ptr_get(seqnum),
                        this->short_header_length_get(seqnum));
      } else
      {
        std::cout.write((const char *) this->header_ptr_get(seqnum),
                        this->header_length_get(seqnum));
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
      assert(currentlength <= sequences_maximum_length);
      if (currentlength == sequences_maximum_length)
      {
        found_maximum_seq_length = true;
      }
#endif
    }
    assert(found_maximum_seq_length);
  }

  /* Overload access operator[], returns
     tuple (header_ptr , header_length, sequence_ptr, sequence_length) */
  std::tuple<const char *, size_t, const char *, size_t> operator[](
      size_t seqnum) const noexcept
  {
    assert(store && seqnum < sequences_number_get());
    const char *head_ptr = static_cast<const char *>(header_ptr[seqnum]);
    const char *seq_ptr = static_cast<const char *>(sequence_ptr[seqnum]);
    const size_t header_len = header_length_get(seqnum);
    const size_t seq_len = sequence_length_get(seqnum);
    return std::tuple<const char *, size_t, const char *, size_t>
                     (head_ptr, header_len, seq_ptr, seq_len);
  }
  Iterator begin()
  {
    return Iterator(sequence_ptr,0);
  }
  Iterator end()
  {
    return Iterator(nullptr,sequences_number);
  }

  template<class T,void (*transformation)(T &,char *,size_t)>
  void transformer(T &t)
  {
    for (size_t snum = 0; snum < sequences_number_get(); snum++)
    {
      transformation(t,sequence_ptr_get(snum),sequence_length_get(snum));
    }
  }
};
#endif  // MULTISEQ_HPP

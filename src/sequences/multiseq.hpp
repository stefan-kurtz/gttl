#ifndef MULTISEQ_HPP
#define MULTISEQ_HPP

#include <tuple>
#include <iostream>
#include <string>
#include <climits>

#include "utilities/gttl_file_open.hpp"
#include "sequences/gttl_seq_iterator.hpp"

/* A class to store various sequences and their header information.
 the inputfile is read in using fasta_reader.
 Constructor may throw std::string
 - std::range_error */

class GttlMultiseq
{
 private:
    size_t sequences_number,
           sequences_total_length,
           sequences_max_length,
           headers_total_length;
    char **sequence_ptr,
         **header_ptr;

    // static counter for all instances of the class
    static uint8_t static_count;
    unsigned char padding_char;

 public:
    /* Constructor
     Inputfile should be in Fasta format, throws std::string */
    GttlMultiseq(const char *inputfile):
      sequences_number(0),
      sequences_total_length(0),
      sequences_max_length(0),
      headers_total_length(0),
      padding_char(UCHAR_MAX - static_count)
    {
      /* set number of instance and increase overall count
         padding_char 122 = 'z', so stop */
      if (padding_char <= 123)
      {
        throw std::string("non letter padding character exhausted");
      }
      static_count++;
      GttlFpType in_fp = gttl_fp_type_open(inputfile, "rb");
      if (in_fp == nullptr)
      {
        throw std::string(": cannot open file");
      }
      const int buf_size = (((size_t) 1) << 14);
      GttlSeqIterator<buf_size> gttl_si_first_pass(in_fp);
      try
      {
        for (auto &&si : gttl_si_first_pass)
        {
          std::string header = std::get<0>(si);
          std::string sequence = std::get<1>(si);

          sequences_total_length += sequence.size();
          assert(header.size() >= 2 &&
                 header[0] == '>' && header[header.size()-1] == '\n');
          headers_total_length += (header.size()-2);
          sequences_number++;
          sequences_max_length = std::max(sequences_max_length,
                                          sequence.size());
        }
      }
      catch(std::string &msg)
      {
        gttl_fp_type_close(in_fp);
        throw msg;
      }
      gttl_fp_type_reset(in_fp);
      sequence_ptr = static_cast<char **>(malloc((sequences_number+1) *
                                                 sizeof *sequence_ptr));
      header_ptr = static_cast<char **>(malloc((sequences_number+1) *
                                               sizeof *header_ptr));

      /* sequence_ptr[0] is malloced to area of size for all symbols
         (sequences_total_length) and padding sumbols, one per sequence,
      therefore sequences_number, plus one for the padding symbol at the start*/
      const size_t sizeof_mem_for_seqs
        = (sequences_total_length + sequences_number + 1) *
           sizeof *sequence_ptr;
      const int return_value
        = posix_memalign(reinterpret_cast<void **>(&sequence_ptr[0]),
                         32, sizeof_mem_for_seqs);
      if (return_value == ENOMEM)
      {
        free(header_ptr);
        free(sequence_ptr);
        gttl_fp_type_close(in_fp);
        StrFormat msg("cannot allocate %lu bytes for all sequences",
                      sizeof_mem_for_seqs);
        throw msg.str();
      }
      sequence_ptr[0][0] = padding_char;
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
      GttlSeqIterator<buf_size> gttl_si_second_pass(in_fp);
      size_t seqnum = 0;
      for (auto &&si : gttl_si_second_pass)
      {
        std::string header = std::get<0>(si);
        std::string sequence = std::get<1>(si);

        memcpy(sequence_ptr[seqnum], sequence.data(), sequence.size());
        *(sequence_ptr[seqnum] + sequence.size()) = padding_char;
        sequence_ptr[seqnum+1] = sequence_ptr[seqnum] + sequence.size() + 1;
        memcpy(header_ptr[seqnum], header.data() + 1, header.size() - 2);
        header_ptr[seqnum+1] = header_ptr[seqnum] + header.size() - 2;
        seqnum++;
      }

      gttl_fp_type_close(in_fp);
    }

    ~GttlMultiseq()
    {
      sequence_ptr[0]--;
      free(sequence_ptr[0]);
      free(header_ptr[0]);
      free(sequence_ptr);
      free(header_ptr);
    }

    size_t number_of_sequences_get(void) const noexcept
    {
      return sequences_number;
    }

    size_t total_length_get(void) const noexcept
    {
      return sequences_total_length;
    }

    unsigned char padding_char_get(void) const noexcept
    {
      return padding_char;
    }

    /* Give the length of sequence seqnum EXCLUDING padding symbol at the end */
    size_t sequence_length_get(size_t seqnum) const noexcept
    {
      assert(seqnum < sequences_number &&
             sequence_ptr[seqnum+1] > sequence_ptr[seqnum]);
      /* To Check whether there are any problems considering the pointer at
      sequences_number goes out of bound and is used for the length of the last
      sequence */
      return (size_t) (sequence_ptr[seqnum+1] - sequence_ptr[seqnum] - 1);
    }

    /* Returns a pointer to the sequence with number seqnum */
    const char *sequence_ptr_get(size_t seqnum) const noexcept
    {
      assert(seqnum < sequences_number);
      return sequence_ptr[seqnum];
    }

    /* Returns length of header. */
    size_t header_length_get(size_t seqnum) const noexcept
    {
      assert(seqnum < sequences_number &&
             header_ptr[seqnum] < header_ptr[seqnum+1]);
      return (size_t) (header_ptr[seqnum+1] - header_ptr[seqnum]);
    }

    /* Returns length of header,
     * short version from start to first space(excluded). */
    size_t short_header_length_get(size_t seqnum) const noexcept
    {
      assert(seqnum < sequences_number &&
             header_ptr[seqnum] < header_ptr[seqnum+1]);
      const char *itr;
      for (itr = header_ptr[seqnum];
           itr < header_ptr[seqnum+1] && !isspace(*itr); itr++)
           /* Nothing */ ;
      return static_cast<size_t>(itr - header_ptr[seqnum]);
    }

    /* Returns a pointer to the header of sequence number seqnum */
    const char *header_ptr_get(size_t seqnum) const noexcept
    {
      assert(seqnum < sequences_number);
      return header_ptr[seqnum];
    }

    size_t maximum_sequence_length_get(void) const
    {
      return sequences_max_length;
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
      size_t seqnum;
#ifndef NDEBUG
      bool found_maximum_seq_length = false;
#endif
      for (seqnum = 0; seqnum < sequences_number; seqnum++)
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
        assert(currentlength <= sequences_max_length);
        if (currentlength == sequences_max_length)
        {
          found_maximum_seq_length = true;
        }
#endif
      }
      assert(found_maximum_seq_length);
    }

    /* Overload access operator[], returns
       tuple (header_ptr , header_length, sequence_ptr, sequence_length) */
    std::tuple<const char*, size_t, const char*, size_t>
          operator[] (size_t seqnum) const noexcept
    {
      assert (seqnum < sequences_number);
      const char *head_ptr = (const char *) header_ptr[seqnum];
      const char *seq_ptr = (const char *) sequence_ptr[seqnum];
      const size_t header_len = header_length_get(seqnum);
      const size_t seq_len = sequence_length_get(seqnum);
      return std::tuple<const char*, size_t, const char*, size_t>
                 (head_ptr, header_len, seq_ptr, seq_len);
    }
};
uint8_t GttlMultiseq::static_count = 0;
#endif  // MULTISEQ_HPP

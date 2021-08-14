#include "sequences/gttl_seq_iterator.hpp"

class GttlMultiseq
{
  private:
    size_t total_length, num_of_sequences, max_sequence_length;
    int max_seq_num_bits,
        max_seq_length_bits;
  public:
    GttlMultiseq(const char *inputfile) :
      total_length(0),
      num_of_sequences(0),
      max_sequence_length(0)
  {
    GttlFpType in_fp = gttl_fp_type_open(inputfile,"rb");
    if (in_fp == nullptr)
    {
      throw std::string(": cannot open file");
    /* check_err.py checked */
    }
    constexpr const int buf_size = 1 << 14;
    GttlSeqIterator<buf_size> gttl_si(in_fp);
    try /* need this, as the catch needs to close the file pointer
           to prevent a memory leak */
    {
      for (auto &&si : gttl_si)
      {
        auto sequence = std::get<1>(si);
        total_length += sequence.size();
        max_sequence_length = std::max(max_sequence_length,sequence.size());
        num_of_sequences++;
      }
    }
    catch (std::string &msg)
    {
      gttl_fp_type_close(in_fp);
      throw msg;
    }
    gttl_fp_type_close(in_fp);
    assert (num_of_sequences > 0);
    max_seq_num_bits = gt_required_bits(num_of_sequences-1);
    max_seq_length_bits = gt_required_bits(max_sequence_length);
  }
  size_t total_length_get(void) const noexcept
  {
    return total_length;
  }
  size_t max_sequence_length_get(void) const noexcept
  {
    return max_sequence_length;
  }
  size_t num_of_sequences_get(void) const noexcept
  {
    return num_of_sequences;
  }
  int max_seq_num_bits_get() const noexcept
  {
    return max_seq_num_bits;
  }
  int max_seq_length_bits_get() const noexcept
  {
    return max_seq_length_bits;
  }
  int sequence_bits_get() const noexcept
  {
    return max_seq_num_bits + max_seq_length_bits;
  }
};

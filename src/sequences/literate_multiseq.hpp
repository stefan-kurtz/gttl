#ifndef LITERATE_MULTISEQ_HPP
#define LITERATE_MULTISEQ_HPP
#include <cstddef>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <iostream>
#include "sequences/alphabet.hpp"
#include "sequences/gttl_multiseq.hpp"
#include "threading/optional_lock.hpp"

template<const char *char_spec,uint8_t undefined_rank>
class LiterateMultiseq
{
  public:
  GttlMultiseq *multiseq;
  static constexpr const GttlAlphabet<char_spec,undefined_rank> alpha{};
  std::array<size_t,alpha.size()+1> rank_dist{};
  private:
  void update_distribution(const char *char_seq,size_t len)
  {
    for (size_t idx = 0; idx < len; idx++)
    {
      rank_dist[static_cast<int>(alpha.char_to_rank(char_seq[idx]))]++;
    }
  }
  void perform_sequence_encoding(char *char_seq,size_t len)
  {
    uint8_t *const encoding = reinterpret_cast<uint8_t *>(char_seq);
    for (size_t idx = 0; idx < len; idx++)
    {
      encoding[idx] = alpha.char_to_rank(char_seq[idx]);
    }
  }
  public:
  LiterateMultiseq(const GttlMultiseq &_multiseq) :
    multiseq(nullptr)
  {
    rank_dist.fill(0);
    for (size_t snum = 0; snum < _multiseq.sequences_number_get(); snum++)
    {
      update_distribution(_multiseq.sequence_ptr_get(snum),
                          _multiseq.sequence_length_get(snum));
    }
  }
  LiterateMultiseq(GttlMultiseq *_multiseq)
    : multiseq(_multiseq)
  {
    rank_dist.fill(0);
    for (size_t snum = 0; snum < multiseq->sequences_number_get(); snum++)
    {
      update_distribution(multiseq->sequence_ptr_get(snum),
                          multiseq->sequence_length_get(snum));
    }
  }
  void perform_sequence_encoding(void)
  {
    for (size_t snum = 0; snum < multiseq->sequences_number_get(); snum++)
    {
      perform_sequence_encoding(multiseq->sequence_ptr_writable_get(snum),
                                multiseq->sequence_length_get(snum));
    }
  }
  void show_rank_dist(std::mutex *cout_mutex) const noexcept
  {
    const OptionalLock lock(cout_mutex);
    for (size_t idx = 0; idx < rank_dist.size(); idx++)
    {
      if (rank_dist[idx] > 0)
      {
        std::cout << "# occurrence\t" << idx << "\t" << rank_dist[idx]
                  << std::endl;
      }
    }
  }
  [[nodiscard]] size_t *rank_dist_copy(void) const noexcept
  {
    size_t *const rdc = reinterpret_cast<size_t *>(
                                 std::malloc(sizeof *rdc * rank_dist.size()));
    std::memcpy(rdc, rank_dist.data(), rank_dist.size() * sizeof *rdc);
    return rdc;
  }
  [[nodiscard]] const std::array<size_t, alpha.size() + 1> &
  rank_dist_get(void) const
  {
    return rank_dist;
  }
  [[nodiscard]] size_t number_of_ranks(void) const noexcept
  {
    return rank_dist.size();
  }
};
#endif

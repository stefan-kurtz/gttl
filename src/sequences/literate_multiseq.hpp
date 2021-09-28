#ifndef LITERATE_MULTISEQ_HPP
#define LITERATE_MULTISEQ_HPP
#include <cstddef>
#include <array>
#include <iostream>
#include "sequences/alphabet.hpp"
#include "sequences/gttl_multiseq.hpp"

template<const char *char_spec,uint8_t undefined_rank>
class LiterateMultiseq
{
  public:
  const GttlMultiseq &multiseq;
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
  public:
  LiterateMultiseq(const GttlMultiseq &_multiseq) :
    multiseq(_multiseq)
  {
    rank_dist.fill(0);
    for (size_t snum = 0; snum < multiseq.sequences_number_get(); snum++)
    {
      update_distribution(multiseq.sequence_ptr_get(snum),
                          multiseq.sequence_length_get(snum));
    }
  }
  void show_rank_dist(void)
  {
    for (size_t idx = 0; idx < rank_dist.size(); idx++)
    {
      if (rank_dist[idx] > 0)
      {
        std::cout << "# occurrence\t" << idx << "\t" << rank_dist[idx]
                  << std::endl;
      }
    }
  }
};
#endif

#ifndef LOC_ALIGN_COORDS_HPP
#define LOC_ALIGN_COORDS_HPP

#include <cstddef>
#include <cstdint>
#include "sequences/gttl_multiseq.hpp"

struct LocalAlignmentCoordinates
{
  size_t ustart, /* start position of local opt. alignment in <useq> */
         usubstringlength, /* length of substring of local opt. align.
                              in <useq> */
         vstart, /* start position of local opt. alignment in <vseq> */
         vsubstringlength; /* length of substring of local opt. align.
                              in <vseq> */
  uint32_t raw_score; /* score of local alignment */
  bool forward_strand;
  LocalAlignmentCoordinates(void) {};
  void show(FILE *fpout,bool dna_alphabet,const GttlMultiseq *multiseq,
            size_t seqnum) const noexcept
  {
    if (usubstringlength + vsubstringlength == 0)
    {
      fprintf(fpout,"%lu\t%lu\t%u",
              ustart,
              vstart,
              raw_score);
    } else
    {
      const size_t vstart_out
        = forward_strand ? vstart
                         : (multiseq->sequence_length_get(seqnum)
                              - vstart
                              - vsubstringlength);
      fprintf(fpout,"%lu\t%lu\t%lu\t%lu\t%u",
              ustart,
              usubstringlength,
              vstart_out,
              vsubstringlength,
              raw_score);
    }
    if (dna_alphabet)
    {
      fprintf(fpout,"\t%c",forward_strand ? '+' : '-');
    }
  }
  bool operator > (const LocalAlignmentCoordinates &other) const noexcept
  {
    return raw_score > other.raw_score or
           (usubstringlength + vsubstringlength > 0 and
            raw_score == other.raw_score and
            usubstringlength + vsubstringlength >
            other.usubstringlength + other.vsubstringlength);
  }
};
#endif

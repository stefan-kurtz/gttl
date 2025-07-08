#ifndef LOC_ALIGN_COORDS_HPP
#define LOC_ALIGN_COORDS_HPP

#include <sys/types.h>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include "sequences/gttl_multiseq.hpp"

struct LocalAlignmentCoordinates
{
  ssize_t ustart = -1; /* start position of local opt. alignment in <useq> */
  ssize_t vstart = -1; /* start position of local opt. alignment in <vseq> */

  size_t usubstringlength = 0, /* length of substring of local opt. align.
                              in <useq> */
         vsubstringlength = 0; /* length of substring of local opt. align.
                              in <vseq> */
  uint32_t raw_score = 0; /* score of local alignment */
  bool forward_strand = true;
  LocalAlignmentCoordinates(void) {};
  void show(FILE *fpout,bool dna_alphabet,const GttlMultiseq *multiseq,
            size_t seqnum) const noexcept
  {
    if (usubstringlength + vsubstringlength == 0)
    {
      fprintf(fpout,"%zu\t%zu\t%u",
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
      fprintf(fpout,"%zu\t%zu\t%zu\t%zu\t%u",
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

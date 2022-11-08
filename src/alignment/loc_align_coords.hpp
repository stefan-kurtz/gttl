#ifndef LOC_ALIGN_COORDS_HPP
#define LOC_ALIGN_COORDS_HPP

#include <cstddef>
#include <cstdint>

struct LocalAlignmentCoordinates
{
  size_t ustart, /* start position of local opt. alignment in <useq> */
         usubstringlength, /* length of substring of local opt. align.
                              in <useq> */
         vstart, /* start position of local opt. alignment in <vseq> */
         vsubstringlength; /* length of substring of local opt. align.
                              in <vseq> */
  uint32_t raw_score; /* score of local alignment */
  void show(FILE *fpout) const noexcept
  {
    fprintf(fpout,"%lu\t%lu\t%lu\t%lu\t%u",
            ustart,
            usubstringlength,
            vstart,
            vsubstringlength,raw_score);
  }
};
#endif

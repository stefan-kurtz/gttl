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
  bool forward_strand;
  void show(FILE *fpout,bool dna_alphabet,const GttlMultiseq *multiseq,
            size_t seqnum) const noexcept
  {
    size_t vstart_out;
    if (forward_strand)
    {
      vstart_out = vstart;
    } else
    {
      const size_t vlen
        = forward_strand ? 0
                         : multiseq->sequence_length_get(seqnum);
      assert(vlen >= vstart + vsubstringlength);
      vstart_out = vlen - vstart - vsubstringlength;
    }
    fprintf(fpout,"%lu\t%lu\t%lu\t%lu\t%u",
            ustart,
            usubstringlength,
            vstart_out,
            vsubstringlength,
            raw_score);
    if (dna_alphabet)
    {
      fprintf(fpout,"\t%c",forward_strand ? '+' : '-');
    }
  }
  bool better(const LocalAlignmentCoordinates &other,
              bool compute_only_end) const noexcept
  {
    return raw_score > other.raw_score ||
           (not compute_only_end &&
            raw_score == other.raw_score &&
            usubstringlength + vsubstringlength >
            other.usubstringlength + other.vsubstringlength);
  }
};
#endif

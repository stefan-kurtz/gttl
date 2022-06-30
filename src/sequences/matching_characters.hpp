#ifndef MATCHING_CHARACTERS_HPP
#define MATCHING_CHARACTERS_HPP
#include "sequences/alphabet.hpp"
#include <cstdbool>

static inline bool matching_characters_wc(char a,char b)
{
  /* This function is only used for self comparison, so we have
     to exclude matches involving wildcards */
  /* if ACGT is considered identical to acgt (no softmasking)
     the low character replace the following by UL_4 */
#ifdef NO_MASKING_LOWERCASE
  static constexpr const alphabet::GttlAlphabet_UL_4 dna_alphabet;
#else
  static constexpr const alphabet::GttlAlphabet_U_4 dna_alphabet;
#endif
  return a == b &&
         dna_alphabet.char_to_rank(a) != dna_alphabet.undefined_rank();
}

static inline bool matching_characters(char a,char b)
{
  return a == b;
}

template<bool (*match_method)(char,char)>
static inline size_t count_mismatches(const char *seq0,const char *seq1,
                                      size_t len)
{
  size_t mismatches = 0;
  for (size_t idx = 0; idx < len; idx++)
  {
    mismatches += (!match_method(seq0[idx],seq1[idx]));
  }
  return mismatches;
}

template<bool (*match_method)(char,char)>
static inline size_t lcplen_fwd(const char *seq0,size_t start0,
                                const char *seq1,size_t start1,
                                GTTL_UNUSED size_t len0,
                                GTTL_UNUSED size_t len1,
                                GTTL_UNUSED size_t seqnum0,
                                GTTL_UNUSED size_t seqnum1)
{
  const char *ptr0 = seq0 + start0,
             *ptr1 = seq1 + start1;
  size_t idx;
  for (idx = 0; match_method(ptr0[idx],ptr1[idx]); idx++)
      /* Nothing */ ;
  return idx;
}

template<bool (*match_method)(char,char)>
static inline size_t lcslen_bwd(const char *seq0,size_t start0,
                                const char *seq1,size_t start1,
                                size_t len0,
                                size_t len1,
                                GTTL_UNUSED size_t seqnum0,
                                GTTL_UNUSED size_t seqnum1)
{
  assert(start0 < len0 && start1 < len1);
  const char *ptr0 = seq0 + len0 - 1 - start0,
             *ptr1 = seq1 + len1 - 1 - start1;
  while (match_method(*ptr0,*ptr1))
  {
    ptr0--;
    ptr1--;
  }
  return static_cast<size_t>((seq0 + len0 - 1 - start0) - ptr0);
}
#endif

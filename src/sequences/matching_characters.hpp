#ifndef MATCHING_CHARACTERS_HPP
#define MATCHING_CHARACTERS_HPP
#include <cstdint>
#include <cstddef>
#include <cassert>
#include "sequences/alphabet.hpp"

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

static inline bool matching_ranks(uint8_t a,uint8_t b)
{
  return a == b;
}

template<bool wildcards>
static inline bool matching_characters_template(char a, char b)
{
  if constexpr (wildcards)
  {
    static constexpr const alphabet::GttlAlphabet_UL_4 dna_alpha{};
    const uint8_t a_rank = dna_alpha.char_to_rank(a);
    const uint8_t b_rank = dna_alpha.char_to_rank(b);
    return a_rank == b_rank && a_rank != dna_alpha.undefined_rank();
  }
  else
  {
    return a == b;
  }
}

template<typename CharType,bool (*match_method)(CharType,CharType)>
static inline size_t count_mismatches(const CharType *seq0,
                                      const CharType *seq1,
                                      size_t len)
{
  size_t mismatches = 0;
  for (size_t idx = 0; idx < len; idx++)
  {
    mismatches += (!match_method(seq0[idx],seq1[idx]));
  }
  return mismatches;
}

template<bool (*match_method)(char,char),bool respect_length>
static inline size_t lcplen_fwd(const char *seq0,size_t start0,
                                const char *seq1,size_t start1,
                                [[maybe_unused]] size_t len0,
                                [[maybe_unused]] size_t len1,
                                [[maybe_unused]] size_t seqnum0,
                                [[maybe_unused]] size_t seqnum1)
{
  if constexpr (respect_length)
  {
    size_t i;
    size_t j;
    for (i = start0, j = start1; i < len0 && j < len1; i++, j++)
    {
      if (not match_method(seq0[i],seq1[j]))
      {
        break;
      }
    }
    return i - start0;
  } else
  {
    size_t idx;
    const char *const ptr0 = seq0 + start0;
    const char *const ptr1 = seq1 + start1;
    for (idx = 0; match_method(ptr0[idx],ptr1[idx]); idx++)
        /* Nothing */;
    return idx;
  }
}

template<bool (*match_method)(char,char)>
static inline size_t lcslen_bwd(const char *seq0,size_t start0,
                                const char *seq1,size_t start1,
                                size_t len0,
                                size_t len1,
                                [[maybe_unused]] size_t seqnum0,
                                [[maybe_unused]] size_t seqnum1)
{
  assert(start0 < len0 && start1 < len1);
  const char *ptr0 = seq0 + len0 - 1 - start0;
  const char *ptr1 = seq1 + len1 - 1 - start1;
  while (match_method(*ptr0,*ptr1))
  {
    ptr0--;
    ptr1--;
  }
  return static_cast<size_t>((seq0 + len0 - 1 - start0) - ptr0);
}

/* The following blockwise comparison function requires a padding of
   the sequences with 7 characters to prevent accessing undefined
   memory. */

static inline size_t block_wise_lcplen_fwd(const char *seq0,size_t start0,
                                           const char *seq1,size_t start1,
                                           [[maybe_unused]] size_t len0,
                                           [[maybe_unused]] size_t len1,
                                           [[maybe_unused]] size_t seqnum0,
                                           [[maybe_unused]] size_t seqnum1)
{
  // Fetch pattern/text blocks
  const uint64_t* block_ptr0 = reinterpret_cast<const uint64_t*>(seq0+start0);
  const uint64_t* block_ptr1 = reinterpret_cast<const uint64_t*>(seq1+start1);
  size_t lcplen;
  // Compare 64-bits blocks
  uint64_t cmp = (*block_ptr0) ^ (*block_ptr1);
  for (lcplen = 0; __builtin_expect(!cmp,0); lcplen += sizeof(uint64_t))
  {
    ++block_ptr0;
    ++block_ptr1;
    cmp = (*block_ptr0) ^ (*block_ptr1);
  }
  // Count equal characters
  const int equal_right_bits = __builtin_ctzl(cmp);
  return lcplen + (equal_right_bits >> 3);
}

static inline size_t block_wise_lcslen_bwd(const char *seq0, size_t start0,
                                           const char *seq1, size_t start1,
                                           size_t len0, size_t len1,
                                           [[maybe_unused]] size_t seqnum0,
                                           [[maybe_unused]] size_t seqnum1)
{
  assert(start0 < len0 && start1 < len1);
  const char *const ptr0     = seq0 + len0 - 8 - start0;
  const char *const ptr1     = seq1 + len1 - 8 - start1;
  const uint64_t *block_ptr0 = reinterpret_cast<const uint64_t *>(ptr0);
  const uint64_t *block_ptr1 = reinterpret_cast<const uint64_t *>(ptr1);

  size_t lcplen;
  uint64_t cmp = (*block_ptr0) ^ (*block_ptr1);
  for (lcplen = 0; __builtin_expect(!cmp, 0); lcplen += sizeof(uint64_t))
  {
    block_ptr0--;
    block_ptr1--;
    cmp = (*block_ptr0) ^ (*block_ptr1);
  }
  const int equal_left_bits = __builtin_clzl(cmp);
  return lcplen + (equal_left_bits >> 3);
}
#endif

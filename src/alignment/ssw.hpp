#ifndef SSW_HPP
#define SSW_HPP
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <cstring>
#include <string>
#include <tuple>
#include <climits>
#include <stdexcept>
#include <format>
#include "utilities/gcc_builtin.hpp"
/* the following is used in sw_simd_uint8|16.hpp */
#include "alignment/simd.hpp"

template<typename Basetype, bool do_print = false>
static void print_simd_int([[maybe_unused]] const char *tag,
                           [[maybe_unused]] const simd_int &vec)
{
  if constexpr (do_print)
  {
    printf("%s", tag);
    simd_int vec_stored;
    simdi_storeu(&vec_stored, vec);
    const Basetype *v = reinterpret_cast<const Basetype *>(&vec_stored);
    for (size_t i = 0; i < 256/(sizeof(Basetype) * CHAR_BIT); i++)
    {
      printf("%d/%d ", static_cast<int>(v[i]), static_cast<Basetype>(v[i]));
    }
    printf("\n");
  }
}

class SimdIntVector
{
  simd_int *ptr;
  size_t num_elements;

  public:
  SimdIntVector(size_t _num_elements)
    : ptr(_num_elements == 0 ? nullptr : new simd_int[_num_elements])
    , num_elements(_num_elements)
  { }

  ~SimdIntVector(void)
  {
    delete[] ptr;
  }

  [[nodiscard]] simd_int *get(void) const noexcept
  {
    return ptr;
  }
  void reset(size_t segment_len)
  {
    assert(ptr != nullptr);
    std::memset(ptr, 0, 4 * segment_len * sizeof *ptr);
  }
  [[nodiscard]] size_t size_in_bytes(void) const noexcept
  {
    return sizeof(*ptr) * num_elements + sizeof(num_elements);
  }
#ifndef NDEBUG
  bool operator==(SimdIntVector &other) const noexcept
  {
    if (num_elements != other.num_elements)
    {
      return false;
    }
    for (size_t idx = 0; idx < num_elements; idx++)
    {
      if (not simd_equal(ptr[idx], other.ptr[idx]))
      {
        return false;
      }
    }
    return true;
  }
#endif
};

// calculate the length including padding required to store len of Basetype
template <typename Basetype>
static inline size_t ssw_len2segment_len(size_t len)
{
  static constexpr const size_t simd_size
    = SIMD_VECSIZE_INT * 4 / sizeof(Basetype);
  return (len + simd_size - 1) / simd_size;
}

class SSWresources
{
  SimdIntVector vectors8, vectors16, vectors32;
#ifndef NDEBUG
  size_t maximum_seq_len;
#endif
  public:
  SSWresources(int which, size_t _maximum_seq_len)
    : vectors8(SimdIntVector(which == 8 or which == 24
                               ? 4 * ssw_len2segment_len<uint8_t>
                                                        (_maximum_seq_len)
                               : 0))
    , vectors16(SimdIntVector(which == 16 || which == 24
                                ? 4 * ssw_len2segment_len<uint16_t>
                                                         (_maximum_seq_len)
                                : 0))
    , vectors32(SimdIntVector(4 *
                              ssw_len2segment_len<uint32_t>(_maximum_seq_len)))
#ifndef NDEBUG
    , maximum_seq_len(_maximum_seq_len)
#endif
  { }
#ifndef NDEBUG
  size_t maximum_seq_len_get(void) const noexcept
  {
    return maximum_seq_len;
  }
#endif
  [[nodiscard]] simd_int *vectors8_get(void) const noexcept
  {
    return vectors8.get();
  }
  [[nodiscard]] simd_int *vectors16_get(void) const noexcept
  {
    return vectors16.get();
  }
  [[nodiscard]] simd_int *vectors32_get(void) const noexcept
  {
    return vectors32.get();
  }
  void reset8(size_t segment_len) { vectors8.reset(segment_len); }
  void reset16(size_t segment_len) { vectors16.reset(segment_len); }
  void reset32(size_t segment_len) { vectors32.reset(segment_len); }
};

class SWsimdResult
{
  public:
  size_t on_dbseq; /* alignment ending position on database */
  size_t on_query; /* alignment ending position on query */
  uint32_t opt_loc_alignment_score;
  SWsimdResult(void)
    : on_dbseq(0)
    , on_query(0)
    , opt_loc_alignment_score(0)
  { }

  SWsimdResult(size_t _on_dbseq, size_t _on_query, uint32_t max_value)
    : on_dbseq(_on_dbseq)
    , on_query(_on_query)
    , opt_loc_alignment_score(max_value)
  { }

  [[nodiscard]] std::string to_string(void) const noexcept
  {
    return
    std::format("on_db_seq={}, on_query_seq={}, opt_loc_alignment_score={}",
                on_dbseq,
                on_query,
                opt_loc_alignment_score);
  }
};

#include "sw_simd_uint8.hpp"
#include "sw_simd_uint16.hpp"
#include "sw_simd_uint32.hpp"

static void ssw_check_valid_sequence(const uint8_t *seq,
                                     size_t seq_len,
                                     size_t alphasize)
{
  for (size_t idx = 0; idx < seq_len; idx++)
  {
    const size_t rank = static_cast<size_t>(seq[idx]);
    if (rank >= alphasize)
    {
      throw std::runtime_error(std::format("illegal rank {} >= {} = alphabet "
                                           "size in position {} of sequence of "
                                           "length {}; make sure that the "
                                           "sequence contain characters for "
                                           "which the pairwise scores are "
                                           "defined",
                                           rank,
                                           alphasize,
                                           idx,
                                           seq_len));
    }
  }
}

template <bool check_valid_sequence, typename Basetype,
          size_t (&sequence_index)(size_t, size_t)>
static inline SimdIntVector ssw_seq_profile(const int8_t *score_vector,
                                            /* this is not used for
                                               sizeof(Basetype) == 2 */
                                            [[maybe_unused]]
                                              uint8_t abs_smallest_score,
                                            size_t alphasize,
                                            const uint8_t *seq,
                                            size_t seq_len)
{
  static_assert(sizeof(Basetype) == 1 or
                sizeof(Basetype) == 2 or
                sizeof(Basetype) == 4);
  if constexpr (check_valid_sequence)
  {
    ssw_check_valid_sequence(seq, seq_len, alphasize);
  }
  static constexpr const size_t simd_size
    = SIMD_VECSIZE_INT * 4 / sizeof(Basetype);
  const size_t segment_len = ssw_len2segment_len<Basetype>(seq_len);
  const SimdIntVector vProfile(alphasize * segment_len);
  Basetype *ptr = reinterpret_cast<Basetype *>(vProfile.get());
  const int8_t *const score_vector_end = score_vector + alphasize * alphasize;

  //  loop over each score row
  for (const int8_t *score_row = score_vector;
       GTTL_IS_LIKELY(score_row < score_vector_end); score_row += alphasize)
  {
    // loop over the simd segments
    for (size_t i = 0; i < segment_len; i++)
    {
      size_t seq_pos = i;

      // loop over each element of a single simd vector
      for (size_t seg_num = 0; GTTL_IS_LIKELY(seg_num < simd_size); seg_num++)
      {
        const int8_t this_score
          = seq_pos < seq_len ? score_row[seq[sequence_index(seq_len, seq_pos)]]
                              : 0;
        if constexpr (sizeof(Basetype) == 1)
        {
          // offset score to ensure it is positive
          *ptr++ = static_cast<Basetype>(abs_smallest_score + this_score);
        } else
        {
          *ptr++ = this_score;
        }
        seq_pos += segment_len;
      }
    }
  }
  return vProfile;
}

static inline size_t sequence_index_forward([[maybe_unused]] size_t seq_len,
                                            size_t seq_pos)
{
  return seq_pos;
}

static inline size_t sequence_index_backward(size_t seq_len, size_t seq_pos)
{
  return seq_len - 1 - seq_pos;
}

#undef PRINT_SSW_PROFILE
#ifdef PRINT_SSW_PROFILE
template<typename T>
static void print_ssw_profile(const SimdIntVector &profile)
{
  constexpr const size_t bits = sizeof(T) * CHAR_BIT;
  printf("Starting to print profile for u%zu.\n#!-", bits);
  for (size_t i = 0; i < profile.size_in_bytes(); i++)
  {
    printf("%d ", reinterpret_cast<T *>(profile.get())[i]);
    constexpr const size_t divisor = size_t(128)/bits;
    if (i % divisor == (divisor - 1))
    {
      printf("\n#!-");
    }
  }
  printf("\n");
}

#define PRINT_ALL_SSW_PROFILES\
        print_profile<uint8_t>(profile_uint8);\
        print_profile<uint16_t>(profile_uint16);\
        print_profile<uint32_t>(profile_uint32)
#else
#define PRINT_ALL_SSW_PROFILES /* Nothing */
#endif

class SSWprofile
{
  const size_t alphasize;
  const int8_t *score_vector; /* flattened scorevector in continues memory */
  const uint8_t abs_smallest_score;
  const uint8_t *query;
  const size_t query_len;
  const SimdIntVector profile_uint8;
  const SimdIntVector profile_uint16;
  const SimdIntVector profile_uint32; /* 0: none */

  public:
  SSWprofile(size_t _alphasize,
             const int8_t *_score_vector,
             int8_t _smallest_score,
             const uint8_t *_query,
             size_t _query_len)
    : alphasize(_alphasize)
    , score_vector(_score_vector)
    , abs_smallest_score(static_cast<uint8_t>(std::abs(_smallest_score)))
    , query(_query)
    , query_len(_query_len)
    , profile_uint8(ssw_seq_profile<true, uint8_t, sequence_index_forward>
                                   (score_vector, abs_smallest_score,
                                    alphasize, query, query_len))
    , profile_uint16(ssw_seq_profile<false, uint16_t, sequence_index_forward>
                                    (score_vector, abs_smallest_score,
                                     alphasize, query, query_len))
    , profile_uint32(ssw_seq_profile<false, uint32_t, sequence_index_forward>
                                    (score_vector, abs_smallest_score,
                                     alphasize, query, query_len))
  {
    PRINT_ALL_SSW_PROFILES;
  }
  ~SSWprofile(void)
  { }
  [[nodiscard]] size_t size_in_bytes(void) const noexcept
  {
    return profile_uint8.size_in_bytes() +
           profile_uint16.size_in_bytes() +
           profile_uint32.size_in_bytes();
  }
  template <bool forward_strand>
  std::tuple<uint32_t, size_t, size_t, size_t, size_t>
    align(SSWresources *ssw_resources,
          const uint8_t *original_dbseq,
          size_t dbseq_len,
          uint8_t weight_gapO,
          uint8_t weight_gapE,
          bool compute_only_end)
  {
    const uint8_t undef_expected_score_uint8 = 0;
    static constexpr const bool forward_reading = true;

    ssw_check_valid_sequence(original_dbseq, dbseq_len, alphasize);
    SWsimdResult forward_ec{
      sw_simd_uint8<forward_reading, forward_strand>
                   (original_dbseq,
                    dbseq_len,
                    dbseq_len,
                    query_len,
                    weight_gapO,
                    weight_gapE,
                    profile_uint8.get(),
                    undef_expected_score_uint8,
                    abs_smallest_score,
                    ssw_resources)};
    int used_uint_bytes = 1;
    if (forward_ec.opt_loc_alignment_score == UINT8_MAX)
    {
      const uint16_t undef_expected_score_uint16 = 0;

      forward_ec = sw_simd_uint16<forward_reading, forward_strand>
                                 (original_dbseq,
                                  dbseq_len,
                                  dbseq_len,
                                  query_len,
                                  weight_gapO,
                                  weight_gapE,
                                  profile_uint16.get(),
                                  undef_expected_score_uint16,
                                  abs_smallest_score,
                                  ssw_resources);
      used_uint_bytes = 2;
      if (forward_ec.opt_loc_alignment_score == INT16_MAX)
      {
        const uint16_t undef_expected_score_uint32 = 0;

        forward_ec = sw_simd_uint32<forward_reading, forward_strand>
                                   (original_dbseq,
                                    dbseq_len,
                                    dbseq_len,
                                    query_len,
                                    weight_gapO,
                                    weight_gapE,
                                    profile_uint32.get(),
                                    undef_expected_score_uint32,
                                    abs_smallest_score,
                                    ssw_resources);
        used_uint_bytes = 4;
        assert(forward_ec.opt_loc_alignment_score < INT32_MAX);
      }
    }
    /* Find the beginning position of the best alignment. */
    size_t ustart;
    size_t usubstringlength;
    size_t vstart;
    size_t vsubstringlength;
    SWsimdResult reverse_ec{};
    if (compute_only_end)
    {
      ustart = forward_ec.on_query;
      usubstringlength = 0;
      vstart = forward_ec.on_dbseq;
      vsubstringlength = 0;
    } else
    {
      if (used_uint_bytes == 1)
      {
        const SimdIntVector vP
          = ssw_seq_profile<false, uint8_t, sequence_index_backward>
                           (score_vector,
                            abs_smallest_score,
                            alphasize,
                            query,
                            forward_ec.on_query + 1);
        reverse_ec = sw_simd_uint8<not forward_reading, forward_strand>
                                  (original_dbseq,
                                   dbseq_len,
                                   forward_ec.on_dbseq + 1,
                                   forward_ec.on_query + 1,
                                   weight_gapO,
                                   weight_gapE,
                                   vP.get(),
                                   static_cast<uint8_t>
                                              (forward_ec
                                                 .opt_loc_alignment_score),
                                   abs_smallest_score,
                                   ssw_resources);
      } else
      {
        if (used_uint_bytes == 2)
        {
          const SimdIntVector vP
            = ssw_seq_profile<false, uint16_t, sequence_index_backward>
                             (score_vector,
                              abs_smallest_score,
                              alphasize,
                              query,
                              forward_ec.on_query + 1);
          reverse_ec = sw_simd_uint16<not forward_reading, forward_strand>
                                     (original_dbseq,
                                      dbseq_len,
                                      forward_ec.on_dbseq + 1,
                                      forward_ec.on_query + 1,
                                      weight_gapO,
                                      weight_gapE,
                                      vP.get(),
                                      forward_ec.opt_loc_alignment_score,
                                      abs_smallest_score,
                                      ssw_resources);
          } else
          {
            assert(used_uint_bytes == 4);
            const SimdIntVector vP
              = ssw_seq_profile<false, uint32_t, sequence_index_backward>
                               (score_vector,
                                abs_smallest_score,
                                alphasize,
                                query,
                                forward_ec.on_query + 1);
            reverse_ec = sw_simd_uint32<not forward_reading, forward_strand>
                                       (original_dbseq,
                                        dbseq_len,
                                        forward_ec.on_dbseq + 1,
                                        forward_ec.on_query + 1,
                                        weight_gapO,
                                        weight_gapE,
                                        vP.get(),
                                        forward_ec.opt_loc_alignment_score,
                                        abs_smallest_score,
                                        ssw_resources);
        }
      }
      assert(forward_ec.on_query >= reverse_ec.on_query);
      ustart = forward_ec.on_query - reverse_ec.on_query;
      usubstringlength = 1 + reverse_ec.on_query;
      vstart = reverse_ec.on_dbseq;
      assert(forward_ec.on_dbseq >= reverse_ec.on_dbseq);
      vsubstringlength = 1 + forward_ec.on_dbseq - reverse_ec.on_dbseq;
    }
    return std::make_tuple(static_cast<uint32_t>
                                      (forward_ec.opt_loc_alignment_score),
                           ustart,
                           usubstringlength,
                           vstart,
                           vsubstringlength);
  }
};
#endif

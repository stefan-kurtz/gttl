/* provide simd method, abstracting from the specific technology
   originally developed by Martin Steinegger, ported to C by
   Leo Foerster and Stefan Kurtz, ported to Neon ny Henning Lindemann. */

#ifndef SIMD_HPP
#define SIMD_HPP
#include <stdlib.h>
#include <stdint.h>

#define AVX512_ALIGN_DOUBLE 64
#define AVX512_VECSIZE_DOUBLE 8
#define AVX512_ALIGN_FLOAT 64
#define AVX512_VECSIZE_FLOAT 16
#define AVX512_ALIGN_INT 64
#define AVX512_VECSIZE_INT 16

#define AVX_ALIGN_DOUBLE 32
#define AVX_VECSIZE_DOUBLE 4
#define AVX_ALIGN_FLOAT 32
#define AVX_VECSIZE_FLOAT 8
#define AVX2_ALIGN_INT 32
#define AVX2_VECSIZE_INT 8

#define SSE_ALIGN_DOUBLE 16
#define SSE_VECSIZE_DOUBLE 2
#define SSE_ALIGN_FLOAT 16
#define SSE_VECSIZE_FLOAT 4
#define SSE_ALIGN_INT 16
#define SSE_VECSIZE_INT 4

#define NEON_ALIGN_INT 16
#define NEON_VECSIZE_INT 4

#define MAX_ALIGN_DOUBLE AVX512_ALIGN_DOUBLE
#define MAX_VECSIZE_DOUBLE AVX512_VECSIZE_DOUBLE
#define MAX_ALIGN_FLOAT AVX512_ALIGN_FLOAT
#define MAX_VECSIZE_FLOAT AVX512_VECSIZE_FLOAT
#define MAX_ALIGN_INT AVX512_ALIGN_INT
#define MAX_VECSIZE_INT AVX512_VECSIZE_INT

#ifdef __AVX2__
#define AVX2
#else

#ifdef __AVX__
#define AVX
#else

#ifdef __SSE4_1__
#define SSE
#else

#ifdef __SSE3__
#define SSE
#else

#ifdef __SSE2__
#define SSE
#else

#ifdef __SSE__
#define SSE
#endif
#endif
#endif
#endif
#endif
#endif

#ifdef __ARM_NEON
#define NEON
#endif

#ifdef NEON
#include <arm_neon.h>
typedef poly128_t simd_int;

#define SIMD_VECSIZE_INT NEON_VECSIZE_INT
#define SIMD_ALIGN_INT NEON_ALIGN_INT

#define simdi16_adds(x, y) \
  vreinterpretq_p128_s16(  \
      vqaddq_s16(vreinterpretq_s16_p128(x), vreinterpretq_s16_p128(y)))
#define simdi16_eq(x, y)  \
  vreinterpretq_p128_u16( \
      vceqq_s16(vreinterpretq_s16_p128(x), vreinterpretq_s16_p128(y)))
#define simdi16_gt(x, y)  \
  vreinterpretq_p128_u16( \
      vcgtq_s16(vreinterpretq_s16_p128(x), vreinterpretq_s16_p128(y)))
#define simdi16_hmax(x) vmaxvq_u16(vreinterpretq_u16_p128(x))
#define simdi16_max(x, y) \
  vreinterpretq_p128_s16( \
      vmaxq_s16(vreinterpretq_s16_p128(x), vreinterpretq_s16_p128(y)))
#define simdi16_set(x) vreinterpretq_p128_s16(vdupq_n_s16(x))
#define simdi32_set(x) vreinterpretq_p128_s32(vdupq_n_s32(x))
#define simdi8_eq(x, y)  \
  vreinterpretq_p128_u8( \
      vceqq_u8(vreinterpretq_u8_p128(x), vreinterpretq_u8_p128(y)))
#define simdi8_hmax(x) vmaxvq_u8(vreinterpretq_u8_p128(x))
#define simdi8_movemask(x) simdi8_movemask_neon_i8(x)
#define simdi8_set(x) vreinterpretq_p128_s8(vdupq_n_s8(x))
#define simdi8_shiftl1(x) \
  vreinterpretq_p128_s8(vextq_s8(vdupq_n_s8(0), vreinterpretq_s8_p128(x), 15))
#define simdi8_shiftl2(x) \
  vreinterpretq_p128_s8(vextq_s8(vdupq_n_s8(0), vreinterpretq_s8_p128(x), 14))
#define simdi_load(x) \
  vldrq_p128((poly128_t*) __builtin_assume_aligned(x, NEON_ALIGN_INT))
#define simdi_loadu(x) vldrq_p128(x)
#define simdi_store(x, y) \
  vstrq_p128((poly128_t*) __builtin_assume_aligned(x, NEON_ALIGN_INT), y)
#define simdi_storeu(x, y) vstrq_p128(x, y)
#define simdui16_subs(x, y) \
  vreinterpretq_p128_u16(   \
      vqsubq_u16(vreinterpretq_u16_p128(x), vreinterpretq_u16_p128(y)))
#define simdui8_adds(x, y) \
  vreinterpretq_p128_u8(   \
      vqaddq_u8(vreinterpretq_u8_p128(x), vreinterpretq_u8_p128(y)))
#define simdui8_max(x, y) \
  vreinterpretq_p128_u8(  \
      vmaxq_u8(vreinterpretq_u8_p128(x), vreinterpretq_u8_p128(y)))
#define simdui8_subs(x, y) \
  vreinterpretq_p128_u8(   \
      vqsubq_u8(vreinterpretq_u8_p128(x), vreinterpretq_u8_p128(y)))

// inspired by
// https://newbedev.com/index.php/
//    sse-mm-movemask-epi8-equivalent-method-for-arm-neon
static int inline simdi8_movemask_neon_i8(poly128_t x)
{
  uint16x8_t temp1 =
      vreinterpretq_u16_u8(vshrq_n_u8(vreinterpretq_u8_p128(x), 7));

  uint32x4_t temp2 = vreinterpretq_u32_u16(vrsraq_n_u16(temp1, temp1, 7));

  uint64x2_t temp3 = vreinterpretq_u64_u32(vrsraq_n_u32(temp2, temp2, 14));

  uint8x16_t temp4 = vreinterpretq_u8_u64(vrsraq_n_u64(temp3, temp3, 28));

  return (int) vgetq_lane_u8(temp4, 0) | ((int) vgetq_lane_u8(temp4, 8)) << 8;
}

#endif

#ifdef AVX512
#include <immintrin.h> /* AVX512 */

/* double support */
#ifndef SIMD_DOUBLE
#define SIMD_DOUBLE
#define ALIGN_DOUBLE AVX512_ALIGN_DOUBLE
#define VECSIZE_DOUBLE AVX512_VECSIZE_DOUBLE
typedef __m512d simd_double;
#define simdf64_add(x, y) _mm512_add_pd(x, y)
#define simdf64_sub(x, y) _mm512_sub_pd(x, y)
#define simdf64_mul(x, y) _mm512_mul_pd(x, y)
#define simdf64_div(x, y) _mm512_div_pd(x, y)
#define simdf64_max(x, y) _mm512_max_pd(x, y)
#define simdf64_load(x) _mm512_load_pd(x)
#define simdf64_store(x, y) _mm512_store_pd(x, y)
#define simdf64_set(x) _mm512_set1_pd(x)
#define simdf64_setzero(x) _mm512_setzero_pd()
#define simdf64_gt(x, y) _mm512_cmpnle_pd_mask(x, y)
#define simdf64_lt(x, y) _mm512_cmplt_pd_mask(x, y)
#define simdf64_or(x, y) _mm512_or_si512(x, y)
#define simdf64_and(x, y) _mm512_and_si512(x, y)
#define simdf64_andnot(x, y) _mm512_andnot_si512(x, y)
#define simdf64_xor(x, y) _mm512_xor_si512(x, y)
#endif

/* float support */
#ifndef SIMD_FLOAT
#define SIMD_FLOAT
#define ALIGN_FLOAT AVX512_ALIGN_FLOAT
#define VECSIZE_FLOAT AVX512_VECSIZE_FLOAT
typedef __m512 simd_float;
#define simdf32_add(x, y) _mm512_add_ps(x, y)
#define simdf32_sub(x, y) _mm512_sub_ps(x, y)
#define simdf32_mul(x, y) _mm512_mul_ps(x, y)
#define simdf32_div(x, y) _mm512_div_ps(x, y)
#define simdf32_rcp(x) _mm512_rcp_ps(x)
#define simdf32_max(x, y) _mm512_max_ps(x, y)
#define simdf32_min(x, y) _mm512_min_ps(x, y)
#define simdf32_load(x) _mm512_load_ps(x)
#define simdf32_store(x, y) _mm512_store_ps(x, y)
#define simdf32_set(x) _mm512_set1_ps(x)
#define simdf32_setzero(x) _mm512_setzero_ps()
#define simdf32_gt(x, y) _mm512_cmpnle_ps_mask(x, y)
#define simdf32_eq(x, y) _mm512_cmpeq_ps_mask(x, y)
#define simdf32_lt(x, y) _mm512_cmplt_ps_mask(x, y)
#define simdf32_or(x, y) _mm512_or_si512(x, y)
#define simdf32_and(x, y) _mm512_and_si512(x, y)
#define simdf32_andnot(x, y) _mm512_andnot_si512(x, y)
#define simdf32_xor(x, y) _mm512_xor_si512(x, y)
#define simdf32_f2i(x) _mm512_cvtps_epi32(x) /* convert float to int */
#define simdf_f2icast(x) _mm512_castps_si512(x)
#endif

/* integer support */
#ifndef SIMD_INT
#define SIMD_INT
#define SIMD_ALIGN_INT AVX512_ALIGN_INT
#define SIMD_VECSIZE_INT AVX512_VECSIZE_INT
typedef __m512i simd_int;
#define simdi32_add(x, y) _mm512_add_epi32(x, y)
#define simdi16_add(x, y) _mm512_add_epi16(x, y)
#define simdi16_adds(x, y) _mm512_adds_epi16(x, y)
#define simdui8_adds(x, y) NOT_YET_IMP()
#define simdi32_sub(x, y) _mm512_sub_epi32(x, y)
#define simdui8_subs(x, y) NOT_YET_IMP()
#define simdi32_mul(x, y) _mm512_mullo_epi32(x, y)
#define simdui8_max(x, y) NOT_YET_IMP()
#define simdi16_max(x, y) _mm512_max_epi32(x, y)
#define simdi32_max(x, y) _mm512_max_epi32(x, y)
#define simdi_load(x) _mm512_load_si512(x)
#define simdi_streamload(x) _mm512_stream_load_si512(x)
#define simdi_store(x, y) _mm512_store_si512(x, y)
#define simdi_storeu(x, y) _mm512_storeu_si512(x, y)
#define simdi32_set(x) _mm512_set1_epi32(x)
#define simdi16_set(x) _mm512_set1_epi16(x)
#define simdi8_set(x) _mm512_set1_epi8(x)
#define simdi_setzero() _mm512_setzero_si512()
#define simdi32_gt(x, y) _mm512_cmpgt_epi32(x, y)
#define simdi8_gt(x, y) NOT_YET_IMP()
#define simdi16_gt(x, y) NOT_YET_IMP()
#define simdi8_eq(x, y) NOT_YET_IMP()
#define simdi32_lt(x, y) _mm512_cmplt_epi32(x, y)
#define simdi16_lt(x, y) NOT_YET_IMP()
#define simdi8_lt(x, y) NOT_YET_IMP()

#define simdi_or(x, y) _mm512_or_si512(x, y)
#define simdi_and(x, y) _mm512_and_si512(x, y)
#define simdi_andnot(x, y) _mm512_andnot_si512(x, y)
#define simdi_xor(x, y) _mm512_xor_si512(x, y)
#define simdi8_shiftl(x, y) NOT_YET_IMP()
#define simdi8_shiftr(x, y) NOT_YET_IMP()
#define simdi8_movemask(x) NOT_YET_IMP()
#define simdi16_extract(x, y) NOT_YET_IMP()
#define simdi16_slli(x, y) _mm512_slli_epi16(x, y) /* shift int left by y */
#define simdi16_srli(x, y) _mm512_srli_epi16(x, y) /* shift int right by y */
#define simdi32_slli(x, y) _mm512_slli_epi32(x, y) /* shift int left by y */
#define simdi32_srli(x, y) _mm512_srli_epi32(x, y) /* shift int right by y */
#define simdi32_i2f(x) _mm512_cvtepi32_ps(x)       /* convert int to float */
#define simdi_i2fcast(x) _mm512_castsi512_ps(x)

#endif /* SIMD_INT */
#endif /* AVX512_SUPPORT */

#ifdef AVX2
/* integer support  (usable with AVX2) */
#ifndef SIMD_INT
#define SIMD_INT
#include <immintrin.h> /* AVX */
#define SIMD_ALIGN_INT AVX2_ALIGN_INT
#define SIMD_VECSIZE_INT AVX2_VECSIZE_INT

inline __m256i _mm256_shift_left1(__m256i a)
{
  __m256i mask = _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0, 0, 3, 0));
  return _mm256_alignr_epi8(a, mask, 16 - 1);
}

inline __m256i _mm256_shift_left2(__m256i a)
{
  __m256i mask = _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0, 0, 3, 0));
  return _mm256_alignr_epi8(a, mask, 16 - 2);
}

inline uint16_t simd_hmax16_avx(const __m256i buffer);
inline uint8_t simd_hmax8_avx(const __m256i buffer);

typedef __m256i simd_int;
#define simdi32_add(x, y) _mm256_add_epi32(x, y)
#define simdi16_add(x, y) _mm256_add_epi16(x, y)
#define simdi16_adds(x, y) _mm256_adds_epi16(x, y)
#define simdui8_adds(x, y) _mm256_adds_epu8(x, y)
#define simdi32_sub(x, y) _mm256_sub_epi32(x, y)
#define simdui16_subs(x, y) _mm256_subs_epu16(x, y)
#define simdui8_subs(x, y) _mm256_subs_epu8(x, y)
#define simdi32_mul(x, y) _mm256_mullo_epi32(x, y)
#define simdi32_max(x, y) _mm256_max_epi32(x, y)
#define simdi16_max(x, y) _mm256_max_epi16(x, y)
#define simdi16_hmax(x) simd_hmax16_avx(x)
#define simdui8_max(x, y) _mm256_max_epu8(x, y)
#define simdi8_hmax(x) simd_hmax8_avx(x)
#define simdi_load(x) _mm256_load_si256(x)
#define simdi_loadu(x) _mm256_loadu_si256(x)
#define simdi_streamload(x) _mm256_stream_load_si256(x)
#define simdi_store(x, y) _mm256_store_si256(x, y)
#define simdi_storeu(x, y) _mm256_storeu_si256(x, y)
#define simdi32_set(x) _mm256_set1_epi32(x)
#define simdi16_set(x) _mm256_set1_epi16(x)
#define simdi8_set(x) _mm256_set1_epi8(x)
#define simdi_setzero() _mm256_setzero_si256()
#define simdi32_gt(x, y) _mm256_cmpgt_epi32(x, y)
#define simdi8_gt(x, y) _mm256_cmpgt_epi8(x, y)
#define simdi16_gt(x, y) _mm256_cmpgt_epi16(x, y)
#define simdi8_eq(x, y) _mm256_cmpeq_epi8(x, y)
#define simdi16_eq(x, y) _mm256_cmpeq_epi16(x, y)
#define simdi32_lt(x, y) _mm256_cmpgt_epi32(y, x) /* inverse */
#define simdi16_lt(x, y) _mm256_cmpgt_epi16(y, x) /* inverse */
#define simdi8_lt(x, y) _mm256_cmpgt_epi8(y, x)
#define simdi_or(x, y) _mm256_or_si256(x, y)
#define simdi_and(x, y) _mm256_and_si256(x, y)
#define simdi_andnot(x, y) _mm256_andnot_si256(x, y)
#define simdi_xor(x, y) _mm256_xor_si256(x, y)
#define simdi8_shiftl1(x) _mm256_shift_left1(x)
#define simdi8_shiftl2(x) _mm256_shift_left2(x)

#define simdi8_shiftr(x, y) _mm256_srli_si256(x, y)
#define simdi8_movemask(x) _mm256_movemask_epi8(x)
#define simdi16_extract(x, y) _mm256_extract_epi16(x, y)
#define simdi16_slli(x, y) _mm256_slli_epi16(x, y) /* shift int left by y */
#define simdi16_srli(x, y) _mm256_srli_epi16(x, y) /* shift int right by y */
#define simdi32_slli(x, y) _mm256_slli_epi32(x, y) /* shift int left by y */
#define simdi32_srli(x, y) _mm256_srli_epi32(x, y) /* shift int right by y */
#define simdi32_i2f(x) _mm256_cvtepi32_ps(x)       /* convert int to float */
#define simdi_i2fcast(x) _mm256_castsi256_ps(x)
#endif /* SIMD_INT */
#endif /* AVX2 */

#ifdef AVX
#include <immintrin.h> /* AVX */
/* double support (usable with AVX1) */
#ifndef SIMD_DOUBLE
#define SIMD_DOUBLE
#define ALIGN_DOUBLE AVX_ALIGN_DOUBLE
#define VECSIZE_DOUBLE AVX_VECSIZE_DOUBLE
typedef __m256d simd_double;
#define simdf64_add(x, y) _mm256_add_pd(x, y)
#define simdf64_sub(x, y) _mm256_sub_pd(x, y)
#define simdf64_mul(x, y) _mm256_mul_pd(x, y)
#define simdf64_div(x, y) _mm256_div_pd(x, y)
#define simdf64_max(x, y) _mm256_max_pd(x, y)
#define simdf64_load(x) _mm256_load_pd(x)
#define simdf64_store(x, y) _mm256_store_pd(x, y)
#define simdf64_set(x) _mm256_set1_pd(x)
#define simdf64_setzero(x) _mm256_setzero_pd()
#define simdf64_gt(x, y) _mm256_cmp_pd(x, y, _CMP_GT_OS)
#define simdf64_lt(x, y) _mm256_cmp_pd(x, y, _CMP_LT_OS)
#define simdf64_or(x, y) _mm256_or_pd(x, y)
#define simdf64_and(x, y) _mm256_and_pd(x, y)
#define simdf64_andnot(x, y) _mm256_andnot_pd(x, y)
#define simdf64_xor(x, y) _mm256_xor_pd(x, y)
#endif /* SIMD_DOUBLE */

/* float support (usable with AVX1) */
#ifndef SIMD_FLOAT
#define SIMD_FLOAT
#define ALIGN_FLOAT AVX_ALIGN_FLOAT
#define VECSIZE_FLOAT AVX_VECSIZE_FLOAT
typedef __m256 simd_float;
#define simdf32_add(x, y) _mm256_add_ps(x, y)
#define simdf32_sub(x, y) _mm256_sub_ps(x, y)
#define simdf32_mul(x, y) _mm256_mul_ps(x, y)
#define simdf32_div(x, y) _mm256_div_ps(x, y)
#define simdf32_rcp(x) _mm256_rcp_ps(x)
#define simdf32_max(x, y) _mm256_max_ps(x, y)
#define simdf32_min(x, y) _mm256_min_ps(x, y)
#define simdf32_load(x) _mm256_load_ps(x)
#define simdf32_store(x, y) _mm256_store_ps(x, y)
#define simdf32_set(x) _mm256_set1_ps(x)
#define simdf32_setzero(x) _mm256_setzero_ps()
#define simdf32_gt(x, y) _mm256_cmp_ps(x, y, _CMP_GT_OS)
#define simdf32_eq(x, y) _mm256_cmp_ps(x, y, _CMP_EQ_OS)
#define simdf32_lt(x, y) _mm256_cmp_ps(x, y, _CMP_LT_OS)
#define simdf32_or(x, y) _mm256_or_ps(x, y)
#define simdf32_and(x, y) _mm256_and_ps(x, y)
#define simdf32_andnot(x, y) _mm256_andnot_ps(x, y)
#define simdf32_xor(x, y) _mm256_xor_ps(x, y)
#define simdf32_f2i(x) _mm256_cvtps_epi32(x)    /* convert float to int */
#define simdf_f2icast(x) _mm256_castps_si256(x) /* compile time cast */
#endif                                          /* SIMD_FLOAT */
#endif                                          /* AVX_SUPPORT */

#ifdef SSE
#include <smmintrin.h> /* SSE4.1 */
/* double support */
#ifndef SIMD_DOUBLE
#define SIMD_DOUBLE
#define ALIGN_DOUBLE SSE_ALIGN_DOUBLE
#define VECSIZE_DOUBLE SSE_VECSIZE_DOUBLE
typedef __m128d simd_double;
#define simdf64_add(x, y) _mm_add_pd(x, y)
#define simdf64_sub(x, y) _mm_sub_pd(x, y)
#define simdf64_mul(x, y) _mm_mul_pd(x, y)
#define simdf64_div(x, y) _mm_div_pd(x, y)
#define simdf64_max(x, y) _mm_max_pd(x, y)
#define simdf64_load(x) _mm_load_pd(x)
#define simdf64_store(x, y) _mm_store_pd(x, y)
#define simdf64_set(x) _mm_set1_pd(x)
#define simdf64_setzero(x) _mm_setzero_pd()
#define simdf64_gt(x, y) _mm_cmpgt_pd(x, y)
#define simdf64_lt(x, y) _mm_cmplt_pd(x, y)
#define simdf64_or(x, y) _mm_or_pd(x, y)
#define simdf64_and(x, y) _mm_and_pd(x, y)
#define simdf64_andnot(x, y) _mm_andnot_pd(x, y)
#define simdf64_xor(x, y) _mm_xor_pd(x, y)
#endif /* SIMD_DOUBLE */

/* float support */
#ifndef SIMD_FLOAT
#define SIMD_FLOAT
#define ALIGN_FLOAT SSE_ALIGN_FLOAT
#define VECSIZE_FLOAT SSE_VECSIZE_FLOAT
typedef __m128 simd_float;
#define simdf32_add(x, y) _mm_add_ps(x, y)
#define simdf32_sub(x, y) _mm_sub_ps(x, y)
#define simdf32_mul(x, y) _mm_mul_ps(x, y)
#define simdf32_div(x, y) _mm_div_ps(x, y)
#define simdf32_rcp(x) _mm_rcp_ps(x)
#define simdf32_max(x, y) _mm_max_ps(x, y)
#define simdf32_min(x, y) _mm_min_ps(x, y)
#define simdf32_load(x) _mm_load_ps(x)
#define simdf32_store(x, y) _mm_store_ps(x, y)
#define simdf32_set(x) _mm_set1_ps(x)
#define simdf32_setzero(x) _mm_setzero_ps()
#define simdf32_gt(x, y) _mm_cmpgt_ps(x, y)
#define simdf32_eq(x, y) _mm_cmpeq_ps(x, y)
#define simdf32_lt(x, y) _mm_cmplt_ps(x, y)
#define simdf32_or(x, y) _mm_or_ps(x, y)
#define simdf32_and(x, y) _mm_and_ps(x, y)
#define simdf32_andnot(x, y) _mm_andnot_ps(x, y)
#define simdf32_xor(x, y) _mm_xor_ps(x, y)
#define simdf32_f2i(x) _mm_cvtps_epi32(x)    /* convert float to int */
#define simdf_f2icast(x) _mm_castps_si128(x) /* compile time cast */
#endif                                       /* SIMD_FLOAT */

/* integer support */
#ifndef SIMD_INT
#define SIMD_INT
#define SIMD_ALIGN_INT SSE_ALIGN_INT
#define SIMD_VECSIZE_INT SSE_VECSIZE_INT
typedef __m128i simd_int;
uint16_t simd_hmax16(const __m128i buffer);
uint8_t simd_hmax8(const __m128i buffer);
#define simdi32_add(x, y) _mm_add_epi32(x, y)
#define simdi16_add(x, y) _mm_add_epi16(x, y)
#define simdi16_adds(x, y) _mm_adds_epi16(x, y)
#define simdui8_adds(x, y) _mm_adds_epu8(x, y)
#define simdi32_sub(x, y) _mm_sub_epi32(x, y)
#define simdui16_subs(x, y) _mm_subs_epu16(x, y)
#define simdui8_subs(x, y) _mm_subs_epu8(x, y)
#define simdi32_mul(x, y) _mm_mullo_epi32(x, y) /* SSE4.1 */
#define simdi32_max(x, y) _mm_max_epi32(x, y)   /* SSE4.1 */
#define simdi16_max(x, y) _mm_max_epi16(x, y)
#define simdi16_hmax(x) simd_hmax16(x)
#define simdui8_max(x, y) _mm_max_epu8(x, y)
#define simdi8_hmax(x) simd_hmax8(x)
#define simdi_load(x) _mm_load_si128(x)
#define simdi_loadu(x) _mm_loadu_si128(x)
#define simdi_streamload(x) _mm_stream_load_si128(x)
#define simdi_storeu(x, y) _mm_storeu_si128(x, y)
#define simdi_store(x, y) _mm_store_si128(x, y)
#define simdi32_set(x) _mm_set1_epi32(x)
#define simdi16_set(x) _mm_set1_epi16(x)
#define simdi8_set(x) _mm_set1_epi8(x)
#define simdi_setzero() _mm_setzero_si128()
#define simdi32_gt(x, y) _mm_cmpgt_epi32(x, y)
#define simdi8_gt(x, y) _mm_cmpgt_epi8(x, y)
#define simdi16_eq(x, y) _mm_cmpeq_epi16(x, y)
#define simdi8_eq(x, y) _mm_cmpeq_epi8(x, y)
#define simdi32_lt(x, y) _mm_cmplt_epi32(x, y)
#define simdi16_lt(x, y) _mm_cmplt_epi16(x, y)
#define simdi8_lt(x, y) _mm_cmplt_epi8(x, y)
#define simdi16_gt(x, y) _mm_cmpgt_epi16(x, y)
#define simdi_or(x, y) _mm_or_si128(x, y)
#define simdi_and(x, y) _mm_and_si128(x, y)
#define simdi_andnot(x, y) _mm_andnot_si128(x, y)
#define simdi_xor(x, y) _mm_xor_si128(x, y)
#define simdi8_shiftl1(x) _mm_slli_si128(x, 1)
#define simdi8_shiftl2(x) _mm_slli_si128(x, 2)
#define simdi8_shiftr(x, y) _mm_srli_si128(x, y)
#define simdi8_movemask(x) _mm_movemask_epi8(x)
#define simdi16_extract(x, y) extract_epi16(x, y)
#define simdi8_extract(x, y) extract_epi8(x, y)
#define simdi16_slli(x, y) _mm_slli_epi16(x, y) /* shift int left by y */
#define simdi16_srli(x, y) _mm_srli_epi16(x, y) /* shift int right by y */
#define simdi32_slli(x, y) _mm_slli_epi32(x, y) /* shift int left by y */
#define simdi32_srli(x, y) _mm_srli_epi32(x, y) /* shift int right by y */
#define simdi32_i2f(x) _mm_cvtepi32_ps(x)       /* convert int to float */
#define simdi_i2fcast(x) _mm_castsi128_ps(x)
#endif /* SIMD_INT */
#endif /* SSE */

#ifdef AVX2
#define SSE_OR_AVX
#endif

#ifdef SSE
#define SSE_OR_AVX
#endif

#ifdef SSE_OR_AVX
inline uint16_t simd_hmax16(const __m128i buffer)
{
  __m128i tmp1 = _mm_subs_epu16(_mm_set1_epi16((short) UINT16_MAX), buffer);
  __m128i tmp3 = _mm_minpos_epu16(tmp1);
  return UINT16_MAX - _mm_cvtsi128_si32(tmp3);
}

inline uint8_t simd_hmax8(const __m128i buffer)
{
  __m128i tmp1 = _mm_subs_epu8(_mm_set1_epi8((char) UINT8_MAX), buffer);
  __m128i tmp2 = _mm_min_epu8(tmp1, _mm_srli_epi16(tmp1, 8));
  __m128i tmp3 = _mm_minpos_epu16(tmp2);
  return (int8_t) (UINT8_MAX - (int8_t) _mm_cvtsi128_si32(tmp3));
}
#endif

#ifdef AVX2
inline uint16_t simd_hmax16_avx(const __m256i buffer)
{
  const __m128i abcd = _mm256_castsi256_si128(buffer);
  const uint16_t first = simd_hmax16(abcd);
  const __m128i efgh = _mm256_extracti128_si256(buffer, 1);
  const uint16_t second = simd_hmax16(efgh);
  return first >= second ? first : second;
}

inline uint8_t simd_hmax8_avx(const __m256i buffer)
{
  const __m128i abcd = _mm256_castsi256_si128(buffer);
  const uint8_t first = simd_hmax8(abcd);
  const __m128i efgh = _mm256_extracti128_si256(buffer, 1);
  const uint8_t second = simd_hmax8(efgh);
  return first >= second ? first : second;
}

/* unclear, why this is implemented by a switch. In case pos <= 15,
   the function is equivalent to _mm256_extract_epi16. Only if pos > 15,
   the function may return 0 which may be different from
   _mm256_extract_epi16. But this could be implemented more
   efficiently by an if then else. */

inline unsigned short extract_epi16(__m256i v, int pos)
{
  switch (pos)
  {
    case 0:
      return _mm256_extract_epi16(v, 0);
    case 1:
      return _mm256_extract_epi16(v, 1);
    case 2:
      return _mm256_extract_epi16(v, 2);
    case 3:
      return _mm256_extract_epi16(v, 3);
    case 4:
      return _mm256_extract_epi16(v, 4);
    case 5:
      return _mm256_extract_epi16(v, 5);
    case 6:
      return _mm256_extract_epi16(v, 6);
    case 7:
      return _mm256_extract_epi16(v, 7);
    case 8:
      return _mm256_extract_epi16(v, 8);
    case 9:
      return _mm256_extract_epi16(v, 9);
    case 10:
      return _mm256_extract_epi16(v, 10);
    case 11:
      return _mm256_extract_epi16(v, 11);
    case 12:
      return _mm256_extract_epi16(v, 12);
    case 13:
      return _mm256_extract_epi16(v, 13);
    case 14:
      return _mm256_extract_epi16(v, 14);
    case 15:
      return _mm256_extract_epi16(v, 15);
  }
  return 0;
}

#else
#ifdef SSE
inline unsigned short extract_epi16(__m128i v, int pos)
{
  switch (pos)
  {
    case 0:
      return _mm_extract_epi16(v, 0);
    case 1:
      return _mm_extract_epi16(v, 1);
    case 2:
      return _mm_extract_epi16(v, 2);
    case 3:
      return _mm_extract_epi16(v, 3);
    case 4:
      return _mm_extract_epi16(v, 4);
    case 5:
      return _mm_extract_epi16(v, 5);
    case 6:
      return _mm_extract_epi16(v, 6);
    case 7:
      return _mm_extract_epi16(v, 7);
  }
  return 0;
}
inline unsigned short extract_epi8(__m128i v, int pos)
{
  switch (pos)
  {
    case 0:
      return _mm_extract_epi8(v, 0);
    case 1:
      return _mm_extract_epi8(v, 1);
    case 2:
      return _mm_extract_epi8(v, 2);
    case 3:
      return _mm_extract_epi8(v, 3);
    case 4:
      return _mm_extract_epi8(v, 4);
    case 5:
      return _mm_extract_epi8(v, 5);
    case 6:
      return _mm_extract_epi8(v, 6);
    case 7:
      return _mm_extract_epi8(v, 7);
    case 8:
      return _mm_extract_epi8(v, 8);
    case 9:
      return _mm_extract_epi8(v, 9);
    case 10:
      return _mm_extract_epi8(v, 10);
    case 11:
      return _mm_extract_epi8(v, 11);
    case 12:
      return _mm_extract_epi8(v, 12);
    case 13:
      return _mm_extract_epi8(v, 13);
    case 14:
      return _mm_extract_epi8(v, 14);
    case 15:
      return _mm_extract_epi8(v, 15);
  }
  return 0;
}
#endif

#endif

static inline bool simd_equal(simd_int a, simd_int b)
{
#ifdef __AVX2__
  const __m256i pcmp = _mm256_cmpeq_epi32(a, b);  // epi8 is fine too
  const unsigned int bitmask = _mm256_movemask_epi8(pcmp);
  return bitmask == 0xffffffffU;
#else
#ifdef NEON
  return a == b;
#else
  static_assert(false);
#endif
#endif
}

#endif /* SIMD_HPP */

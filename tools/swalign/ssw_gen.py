#!/usr/bin/env python3

import sys
import argparse
import jinja2
# https://realpython.com/primer-on-jinja-templating/#install-jinja

environment = jinja2.Environment()
template = environment.from_string('''#include <alignment/ssw.hpp>
template<bool forward_reading,bool forward_strand> static SWsimdResult sw_simd_uint{{ width }} (
                                 GTTL_UNUSED const uint8_t *original_dbseq,
                                 GTTL_UNUSED size_t original_dbseq_len,
                                 size_t dbseq_len,
                                 size_t query_len,
                                 const uint8_t weight_gapO,
                                 const uint8_t weight_gapE,
                                 const simd_int *vProfile,
                                 /* best alignment score: used to terminate
                                    the matrix calculation when locating the
                                    alignment beginning point. If this score
                                    is set to 0, it will not be used */
                                 uint{{ width }}_t expected_score,
                                 uint8_t abs_smallest_score,
                                 SSWresources *ssw_resources)
{
#if SSW_SIMD_DEBUG > 0
  uint64_t column_count = 0;
  uint64_t column_max_move_count = 0;
#endif
  SWsimdResult sw_simd_result(0,query_len - 1,UINT{{ width }}_MAX);
  const size_t simd_size = SIMD_VECSIZE_INT * {{ 32//width }},
               segment_len = (query_len + simd_size - 1) / simd_size;
  uint{{ width }}_t max_align_score = 0;
  int step;
  bool own_resources;
  const simd_int vZero = simdi32_set(0),
                 vGapO = simdi{{ width }}_set(weight_gapO),
                 vGapE = simdi{{ width }}_set(weight_gapE);{{ v_bias_init_var }}
  print_simd_int<uint{{ width }}_t>("vGapO: ", vGapO);
  simd_int vTemp;
  uint32_t cmp;

  if (ssw_resources == NULL)
  {
    ssw_resources = new SSWresources({{ width }}, query_len);
    own_resources = true;
  } else
  {
    assert(query_len <= ssw_resources->maximum_seq_len_get());
    ssw_resources->reset{{ width }}(segment_len);
    own_resources = false;
  }
  simd_int *pvHStore = ssw_resources->vectors{{ width }}_get();
  simd_int *pvHLoad = pvHStore + segment_len;
  simd_int *pvE = pvHLoad + segment_len;
  simd_int *pvHmax = pvE + segment_len;
  simd_int *pvHStoreNext = pvHLoad;
  simd_int *pvHStoreNextNext = pvHStore;

  int64_t dbseq_pos, dbseq_pos_end;
  assert(dbseq_len > 0);
  if constexpr (forward_reading)
  {
    dbseq_pos = 0;
    dbseq_pos_end = static_cast<int64_t>(dbseq_len);
    step = 1;
  } else
  {
    dbseq_pos = static_cast<int64_t>(dbseq_len - 1);
    dbseq_pos_end = -1;
    step = -1;
  }

  /* outer loop to process the database sequence */
  while (dbseq_pos != dbseq_pos_end)
  {
#if SSW_SIMD_DEBUG > 0
    column_count++;
#endif
    assert(dbseq_pos >= 0);
    uint8_t current_char;
    if constexpr (forward_strand)
    {
      current_char = original_dbseq[dbseq_pos];
    } else
    {
      if constexpr (forward_reading)
      {
        assert(dbseq_len == original_dbseq_len);
      }
      current_char = original_dbseq[original_dbseq_len - 1 - dbseq_pos];
      current_char = complement_uint8_wc_remains(current_char);
    }
    //printf(\"Current_char %d\\n\", current_char);
    size_t segment_pos;
    simd_int e,
             *pv,
             vF = vZero,
             vMaxColumn = vZero,
             vH = pvHStore[segment_len - 1];
    const simd_int *vP = vProfile + segment_len * static_cast<size_t>(current_char);

    print_simd_int<uint{{ width }}_t>("Initial vH: ", vH);
    vH = simdi8_shiftl{{width//8}}(vH); /* Shift the value in vH left by 2 byte. */

    print_simd_int<uint{{ width }}_t>("vH shifted: ", vH);

    /* Swap the 2 H buffers. */
    pvHLoad = pvHStore;
    pvHStore = pvHStoreNext;
    pv = pvHStoreNext;
    pvHStoreNext = pvHStoreNextNext;
    pvHStoreNextNext = pv;

    /* inner loop to process the query sequence */
    for (segment_pos = 0; GTTL_IS_LIKELY(segment_pos < segment_len); ++segment_pos)
    {
      vH = simd{{ ui_for_8 }}{{ width }}_adds(vH, simdi_load(vP + segment_pos));{{ vH_subs }}
      print_simd_int<uint{{ width }}_t>("for loop 1 vH: ", vH);

      /* Get max from vH, vE and vF. */
      e = simdi_load(pvE + segment_pos);
      vH = simd{{ ui_for_8 }}{{ width }}_max(vH, e);
      vH = simd{{ ui_for_8 }}{{ width }}_max(vH, vF);
      vMaxColumn = simd{{ ui_for_8 }}{{ width }}_max(vMaxColumn, vH);

      print_simd_int<uint{{ width }}_t>("for loop 2 vH: ", vH);
      print_simd_int<uint{{ width }}_t>("for loop 2 vMaxColumn: ", vMaxColumn);

      /* Save vH values. */
      simdi_store(pvHStore + segment_pos, vH);

      /* Update vE value. */
      vH = simdui{{ width }}_subs(vH, vGapO); /* saturation arithmetic, result >= 0 */
      e = simdui{{ width }}_subs(e, vGapE);
      e = simd{{ ui_for_8 }}{{ width }}_max(e, vH);
      simdi_store(pvE + segment_pos, e);

      /* Update vF value. */{{ correct_comment }}
      vF = simdui{{ width }}_subs(vF, vGapE);
      vF = simd{{ ui_for_8 }}{{ width }}_max(vF, vH);

      /* Load the next vH. */
      vH = simdi_load(pvHLoad + segment_pos);
    }

    /* Lazy_F loop: has been revised to disallow adjacent insertion and
       then deletion, so do not update E(i, segment_pos), learn from SWPS3 */
#ifndef AVX2
#define SSW_MAX_CMP_VALUE UINT16_MAX
#else
#define SSW_MAX_CMP_VALUE UINT32_MAX
#endif /* AVX2 */
#if 8 == {{ width }}

    /* reset pointers to the start of the saved data */
    vH = simdi_load(pvHStore);

    /* the computed vF value is for the given column.  since */
    /* we are at the end, we need to shift the vF value over */
    /* to the next column. */
    vF = simdi8_shiftl1(vF);
    vTemp = simdui8_subs(vH, vGapO);
    vTemp = simdui8_subs(vF, vTemp);
    vTemp = simdi8_eq(vTemp, vZero);

    for (cmp = simdi8_movemask(vTemp), segment_pos = 0; cmp != SSW_MAX_CMP_VALUE;
         cmp = simdi8_movemask(vTemp))
    {
      vH = simdui8_max (vH, vF);
      vMaxColumn = simdui8_max(vMaxColumn, vH);
      simdi_store(pvHStore + segment_pos, vH);
      vF = simdui8_subs (vF, vGapE);
      segment_pos++;
      if (segment_pos >= segment_len)
      {
        segment_pos = 0;
        vF = simdi8_shiftl1(vF);
      }
      vH = simdi_load(pvHStore + segment_pos);
      vTemp = simdui8_subs(vH, vGapO);
      vTemp = simdui8_subs(vF, vTemp);
      vTemp = simdi8_eq(vTemp, vZero);
    }
#else
    for (size_t k = 0; GTTL_IS_LIKELY(k < simd_size); ++k)
    {
      vF = simdi8_shiftl{{ width//8 }}(vF);
      for (segment_pos = 0; GTTL_IS_LIKELY(segment_pos < segment_len); ++segment_pos)
      {
        vH = simdi_load(pvHStore + segment_pos);
        vH = simdi{{ width }}_max(vH, vF);
        vMaxColumn = simdi{{ width }}_max(vMaxColumn, vH); /*newly added line */
        simdi_store(pvHStore + segment_pos, vH);
        /* Update vF value. */
        vH = simdui{{ width }}_subs(vH, vGapO);
        vF = simdui{{ width }}_subs(vF, vGapE);
        if (GTTL_IS_UNLIKELY(!simdi8_movemask(simdi{{ width }}_gt(vF, vH))))
        {
          break;
        }
      }
      if (segment_pos < segment_len)
      {
        break;
      }
    }
#endif

    simd_int vMaxScore = vZero; /* highest score of the whole matrix. */
    print_simd_int<uint{{ width }}_t>("vMaxScore: ", vMaxScore);
    print_simd_int<uint{{ width }}_t>("vMaxColumn: ", vMaxColumn);
    vMaxScore = simd{{ ui_for_8 }}{{ width }}_max(vMaxScore, vMaxColumn);
    print_simd_int<uint{{ width }}_t>("vMaxScore: ", vMaxScore);
    simd_int vMaxMark = vZero;  /* highest score until previous column. */
    vTemp = simdi{{ width }}_eq(vMaxMark, vMaxScore);
    cmp = simdi8_movemask(vTemp);
    if (cmp != SSW_MAX_CMP_VALUE)
    {
      const uint{{ width }}_t local_max_score = simdi{{ width }}_hmax(vMaxScore);

      vMaxMark = vMaxScore;
      if (GTTL_IS_LIKELY(local_max_score > max_align_score))
      {
        max_align_score = local_max_score;
        if (static_cast<uint32_t>(max_align_score) +
            static_cast<uint32_t>(abs_smallest_score) >= UINT{{ width }}_MAX)
        {
          break;  /*overflow */
        }
        sw_simd_result.on_dbseq = static_cast<size_t>(dbseq_pos);
        /* Store the column with the highest alignment score in order to
           trace the alignment ending position on query. */
        // memcpy(pvHmax,pvHStore,segment_len * sizeof *pvHmax);
#if SSW_SIMD_DEBUG > 0
        column_max_move_count++;
#endif
        pv = pvHmax;
        pvHmax = pvHStoreNextNext;
        pvHStoreNextNext = pv;
      }
    }

    /* Record the max score of current column. */
    if (expected_score > 0 && simdi{{ width }}_hmax(vMaxColumn) == expected_score)
    {
      break;
    }
    dbseq_pos += step;
#if SSW_SIMD_DEBUG > 1
    uint{{ width }}_t *ptr = reinterpret_cast<uint{{ width }}_t *>(pvHStore);
    for (size_t i = 0; i < segment_len * simd_size; i++) {
      printf("%4d", ptr[(i % segment_len) * simd_size + i / segment_len]);
    }
    printf("\\n");
#endif
  }

  if (static_cast<uint32_t>(max_align_score) + static_cast<uint32_t>(abs_smallest_score) < UINT{{ width }}_MAX)
  {
    /* Trace the column with the max alignment score for the ending
       position on query. */
    uint{{ width }}_t *ptr = reinterpret_cast<uint{{ width }}_t*>(pvHmax);
    const size_t column_len = segment_len * simd_size;
    for (size_t i = 0; GTTL_IS_LIKELY(i < column_len); ++i, ++ptr)
    {
      if (*ptr == max_align_score)
      {
        const size_t current_end
          = i / simd_size + (i % simd_size) * segment_len;
        if (current_end < sw_simd_result.on_query)
        {
          sw_simd_result.on_query = current_end;
        }
      }
    }
    sw_simd_result.opt_loc_alignment_score = max_align_score;
  }
  if (own_resources)
  {
    delete ssw_resources;
  }
#if SSW_SIMD_DEBUG > 0
  printf("alignment uint{{ width }} %zu/%zu\\n", column_max_move_count, column_count);
#endif
  return sw_simd_result;
}''')


def parse_arguments(argv):
    p = argparse.ArgumentParser(description='generate code for SSW-functions')
    p.add_argument('width', type=int, choices=[8, 16, 32])
    return p.parse_args(argv)


args = parse_arguments(sys.argv[1:])

if args.width == 8:
    v_bias_init_expr = ('\nconst simd_int vBias = simdi{}_set(abs_smallest_score);'
                        .format(args.width))
else:
    v_bias_init_expr = ''
print('/* generated by {} DO NOT EDIT */'.format(' '.join(sys.argv)))
print('#ifndef SW_SIMD_UINT{}_HPP'.format(args.width))
print('#define SW_SIMD_UINT{}_HPP'.format(args.width))
print(template.render(width=args.width,
                      ui_for_8='ui' if args.width == 8 else 'i',
                      vH_subs='\nvH = simdui8_subs(vH, vBias); /* vH will be always > 0 */' if args.width == 8 else '',
                      v_bias_init_var=v_bias_init_expr,
                      correct_comment='\n/* correct 8 -> 16 in the next line */' if args.width >= 16 else ''))
print('#endif')

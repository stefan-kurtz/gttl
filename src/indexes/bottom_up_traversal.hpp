#ifndef BOTTOM_UP_TRAVERSAL_HPP
#define BOTTOM_UP_TRAVERSAL_HPP
#include <cstddef>
#include <cstdint>
#include <cstdbool>
#include <vector>
#include <tuple>
#include <algorithm>
#include <climits>
#include "utilities/bitpacker.hpp"
#include "utilities/bytes_unit.hpp"

/* This file implements the essence of the bottom up traversal, independent
   of the specific needs for the maxpair computation. */

template<int sizeof_unit>
static std::pair<uint32_t,uint32_t> decode_seqnum_relpos(
                                         const BytesUnit<sizeof_unit,2> &suffix,
                                         const GttlBitPacker<sizeof_unit,2> &bp)
{
  const uint64_t leafnumber_seqnum = suffix.template decode_at<0>(bp);
  const uint64_t leafnumber_relpos = suffix.template decode_at<1>(bp);
  assert(leafnumber_relpos <= UINT32_MAX && leafnumber_seqnum <= UINT32_MAX);
  return {static_cast<uint32_t>(leafnumber_seqnum),
          static_cast<uint32_t>(leafnumber_relpos)};
}

template<typename Basetype>
class BottomUpTraversalStack
{
  private:
  Basetype *space;
  size_t allocated, nextfree;
  public:
  BottomUpTraversalStack(void) :
    space(nullptr),
    allocated(0),
    nextfree(0) {}
  ~BottomUpTraversalStack(void)
  {
    free(space);
  }
  void push_back(size_t lcp, size_t lb)
  {
    if (nextfree >= allocated)
    {
      size_t previously_allocated = allocated;
      allocated = (allocated * 1.2) + 32;
      space = static_cast<Basetype *>(realloc(space,allocated * sizeof *space));
      for (size_t idx = previously_allocated; idx < allocated; idx++)
      {
        space[idx].lcp = 0;
        space[idx].lb = 0;
        space[idx].rb = ULONG_MAX;
      }
    }
    space[nextfree].lcp = lcp;
    space[nextfree].lb = lb;
    nextfree++;
  }
  Basetype *back_ptr(void) const noexcept
  {
    assert(nextfree > 0);
    return space + nextfree - 1;
  }
  void pop_back(void)
  {
    assert(nextfree > 0);
    nextfree--;
  }
  size_t size(void) const noexcept
  {
    return nextfree;
  }
};

template<typename StateType,typename IntervalRecord>
using ProcessLeafEdgeFunction = void (*)(StateType *,
                                         bool,
                                         size_t,
                                         IntervalRecord *,
                                         uint32_t,
                                         uint32_t,
                                         bool);

template<typename StateType,typename IntervalRecord>
using ProcessBranchingEdgeFunction = void (*)(StateType *,
                                              bool,
                                              size_t,
                                              size_t,
                                              IntervalRecord *,
                                              size_t,
                                              size_t,
                                              size_t,
                                              IntervalRecord *,
                                              bool);

template<typename Basetype>
struct BUItvinfo
{
  size_t lcp, lb, rb;
  Basetype info;
};

template<class SuffixArray,
         typename StateType,typename IntervalRecord,int sizeof_unit>
static void bottomup_generic(bool with_mmap,
                             StateType *bu_state,
                             const SuffixArray *suffixarray,
                             size_t nonspecial_suffixes,
                             ProcessLeafEdgeFunction<StateType,IntervalRecord>
                               process_leafedge,
                             ProcessBranchingEdgeFunction<StateType,
                                                          IntervalRecord>
                               process_branchingedge)
{
  const int sequences_number_bits = suffixarray->sequences_number_bits_get();
  const int sequences_length_bits = suffixarray->sequences_length_bits_get();
  int first_group_bits;
  bool first_edge_from_root = true;
  BUItvinfo<IntervalRecord> *last_interval = nullptr;
  BottomUpTraversalStack<BUItvinfo<IntervalRecord>> stack{};

  assert(sequences_number_bits + sequences_length_bits <=
         sizeof_unit * CHAR_BIT);
  if (suffixarray->sequences_number_get() == 1)
  {
    assert(sequences_number_bits == 0);
    first_group_bits = sizeof_unit * CHAR_BIT - sequences_length_bits;
  } else
  {
    assert(sequences_number_bits > 0);
    first_group_bits = sequences_number_bits;
  }
  const GttlBitPacker<sizeof_unit,2> bp({first_group_bits,
                                         sequences_length_bits});
  const BytesUnit<sizeof_unit,2> *bu_suftab;
  if (with_mmap)
  {
    const uint8_t* suftab_bytes = suffixarray->get_mmap_suftab_bytes();
    bu_suftab = reinterpret_cast<const BytesUnit<sizeof_unit,2> *>
                                (suftab_bytes);
  } else
  {
    const std::vector<uint8_t> &suftab_bytes = suffixarray->get_suftab_bytes();
    bu_suftab = reinterpret_cast<const BytesUnit<sizeof_unit,2> *>
                                (suftab_bytes.data());
  }
  const LCPtable &lcptable = suffixarray->lcptable_get();
  stack.push_back(0,0);
  size_t interval_bound = 0;
  for (auto it = lcptable.begin(); it != lcptable.end(); ++it)
  {
    const unsigned int lcpvalue = *it; /* at interval_bound + 1 */
    auto seqnum_relpos = decode_seqnum_relpos(bu_suftab[interval_bound],bp);
    assert(stack.size() > 0);
    if (lcpvalue <= stack.back_ptr()->lcp)
    {
      bool first_edge;
      const bool last_child
        = static_cast<bool>(lcpvalue < stack.back_ptr()->lcp);
      if (stack.back_ptr()->lcp > 0 || !first_edge_from_root)
      {
        first_edge = false;
      } else
      {
        first_edge = true;
        first_edge_from_root = false;
      }
      process_leafedge(bu_state,
                       first_edge,
                       stack.back_ptr()->lcp,
                       &stack.back_ptr()->info,
                       std::get<0>(seqnum_relpos),
                       std::get<1>(seqnum_relpos),
                       last_child);
    }
    assert(last_interval == nullptr);
    while (lcpvalue < stack.back_ptr()->lcp)
    {
      last_interval = stack.back_ptr();
      stack.pop_back();
      last_interval->rb = interval_bound;
      if (lcpvalue <= stack.back_ptr()->lcp)
      {
        bool first_edge;
        if (stack.back_ptr()->lcp > 0 || !first_edge_from_root)
        {
          first_edge = false;
        } else
        {
          first_edge = true;
          first_edge_from_root = false;
        }
        const bool last_child
          = static_cast<bool>(lcpvalue < stack.back_ptr()->lcp);
        process_branchingedge(bu_state,
                              first_edge,
                              stack.back_ptr()->lcp,
                              stack.back_ptr()->lb,
                              &stack.back_ptr()->info,
                              last_interval->lcp,
                              last_interval->lb,
                              last_interval->rb,
                              &last_interval->info,
                              last_child);
        last_interval = nullptr;
      }
    }
    if (lcpvalue > stack.back_ptr()->lcp)
    {
      if (last_interval != nullptr)
      {
        const size_t last_interval_lb = last_interval->lb,
                     last_interval_lcp = last_interval->lcp,
                     last_interval_rb = last_interval->rb;
        stack.push_back(lcpvalue,last_interval_lb);
        static constexpr const bool first_edge = true;
        process_branchingedge(bu_state,
                              first_edge,
                              stack.back_ptr()->lcp,
                              stack.back_ptr()->lb,
                              &stack.back_ptr()->info,
                              last_interval_lcp,
                              last_interval_lb,
                              last_interval_rb,
                              nullptr, /* child not used for first_edge=true */
                              false);
        last_interval = nullptr;
      } else
      {
        stack.push_back(lcpvalue,interval_bound);
        static constexpr const bool first_edge = true;
        process_leafedge(bu_state,
                         first_edge,
                         stack.back_ptr()->lcp,
                         &stack.back_ptr()->info,
                         std::get<0>(seqnum_relpos),
                         std::get<1>(seqnum_relpos),
                         false);
      }
    }
    interval_bound++;
    if (interval_bound == nonspecial_suffixes)
    {
      break;
    }
  }
  assert(stack.size() > 0 && stack.back_ptr()->lcp == 0);
}
#endif

#ifndef BOTTOM_UP_TRAVERSAL_HPP
#define BOTTOM_UP_TRAVERSAL_HPP
#include <cstddef>
#include <cstdint>
#include <cstdbool>
#include <climits>

/* This file implements the essence of the bottom up traversal, independent
   of the specific needs for the maxpair computation. */

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

template<class StateClass,class IntervalRecord,typename SuftabType>
using ProcessLeafEdgeFunction = void (*)(StateClass *,
                                         bool,
                                         size_t,
                                         IntervalRecord *,
                                         const SuftabType &,
                                         bool);

template<class StateClass,class IntervalRecord>
using ProcessBranchingEdgeFunction = void (*)(StateClass *,
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

template<class SuftabClass, class LCPtabClass,
         class StateClass,class IntervalRecord, typename SuftabType>
static void bottomup_generic(StateClass *bu_state,
                             const SuftabClass &suftab,
                             const LCPtabClass &lcptab,
                             ProcessLeafEdgeFunction<StateClass,IntervalRecord,
                                                     SuftabType>
                               process_leafedge,
                             ProcessBranchingEdgeFunction<StateClass,
                                                          IntervalRecord>
                               process_branchingedge)
{
  bool first_edge_from_root = true;
  BUItvinfo<IntervalRecord> *last_interval = nullptr;
  BottomUpTraversalStack<BUItvinfo<IntervalRecord>> stack{};

  stack.push_back(0,0);
  size_t interval_bound = 0;
  const size_t nonspecial_suffixes = suftab.nonspecial_suffixes_get();
  for (auto it = lcptab.begin(); it != lcptab.end(); ++it)
  {
    const unsigned int lcpvalue = *it; /* at interval_bound + 1 */
    auto seqnum_relpos = suftab[interval_bound];
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
                       seqnum_relpos,
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
                         seqnum_relpos,
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

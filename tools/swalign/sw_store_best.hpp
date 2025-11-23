#ifndef SW_STORE_BEST_HPP
#define SW_STORE_BEST_HPP
#include <cstddef>
#include <cassert>
#include <vector>
#include "utilities/fs_prio_store.hpp"
#include "alignment/loc_align_coords.hpp"

struct StoredLocalAlignmentCoordinates
{
  LocalAlignmentCoordinates coords;
  size_t u_seqnum, v_seqnum;
  StoredLocalAlignmentCoordinates(void)
    : coords({})
    , u_seqnum(0)
    , v_seqnum(0)
  {}
  StoredLocalAlignmentCoordinates(const LocalAlignmentCoordinates &_coords,
                                  size_t _u_seqnum,
                                  size_t _v_seqnum)
    : coords(_coords)
    , u_seqnum(_u_seqnum)
    , v_seqnum(_v_seqnum)
  {}
  bool operator > (const StoredLocalAlignmentCoordinates &other) const
  {
    return coords > other.coords;
  }
};

using SWResultVector = FSPrioStore<StoredLocalAlignmentCoordinates>;

class SWStoreBestResultsShared
{
  public:
  SWStoreBestResultsShared(void) = default;
  bool process(SWResultVector *store,
               const LocalAlignmentCoordinates &best_coords,
               size_t i, size_t j) const
  {
    store->add(StoredLocalAlignmentCoordinates(best_coords,i,j));
    return false;
  }
};

class SWStoreBestResultsGetThreadRelated
{
  std::vector<SWResultVector *> results;
  public:
  SWStoreBestResultsGetThreadRelated(size_t num_threads,size_t best)
    : results({})
  {
    for (size_t t = 0; t < num_threads; t++)
    {
      results.push_back(new SWResultVector(best));
    }
  }
  ~SWStoreBestResultsGetThreadRelated(void)
  {
    for (auto v : results)
    {
      delete v;
    }
  }
  SWResultVector *operator [](size_t t) const
  {
    assert(t < results.size());
    return results[t];
  }
};
#endif

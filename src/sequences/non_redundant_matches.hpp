#ifndef NON_REDUNDANT_MATCHES_HPP
#define NON_REDUNDANT_MATCHES_HPP
#include <cstdint>
#include <utility>
#include <cstddef>
#include <vector>
#include <set>
#include <algorithm>
#include <cassert>
#include "utilities/multibitvector.hpp"

template <typename ElementClass>
class NonRedundantMatches
{
  class CmpIds
  {
    const std::vector<ElementClass> &elements;

    public:
    CmpIds(const std::vector<ElementClass> &_elements)
      : elements(_elements)
    {}
    // comparison used to define order in BST
    bool operator ()(const uint32_t i, const uint32_t j) const noexcept
    {
      return elements[i].superior_weight_tie_primary_startpos(elements[j]);
    }
  };

  static constexpr const bool track_count = false;
  Multibitvector<track_count> good;
  // 1-bit for 'good' elements and 0-bit for 'bad' elements

  public:
  NonRedundantMatches(const std::vector<ElementClass> &elements)
    : good(elements.size())
  {
    if (elements.size() <= 1)
    {
      if (elements.size() == 1)
      {
        good.set(0);
      }
      return;
    }
    // BST of segment ids, ordered by their weights
    const CmpIds cmp_ids(elements);
    std::set<uint32_t, CmpIds> status(cmp_ids);

    // vector of pairs <endposition, id incl. event type>
    // id is *even* for 'sequence start' and *odd* for 'sequence end'
    using Event = std::pair<uint32_t,uint32_t>;
    std::vector<Event> event_schedule{};
    event_schedule.reserve(2 * elements.size());
    assert(((elements.size()-1) << 1) + 1 <= UINT32_MAX);
    for (uint32_t idx = 0; idx < static_cast<uint32_t>(elements.size()); idx++)
    {
      const auto primary_startpos = elements[idx].primary_startpos_get();
      event_schedule.emplace_back(primary_startpos,idx << 1);
      event_schedule.emplace_back(primary_startpos +
                                  elements[idx].primary_len_get(),
                                  (idx << 1) + 1);
    }
    std::sort(event_schedule.begin(), event_schedule.end());

    for (auto it_sweep_line = event_schedule.begin();
         it_sweep_line != event_schedule.end(); /* Nothing */)
    {
      auto current = std::get<0>(*it_sweep_line);
      auto it_differ = std::find_if(it_sweep_line, event_schedule.end(),
                                    [&](const Event &e)
                                    {
                                      return std::get<0>(e) != current;
                                    });
      /*
      idx_differ = idx_sweep_line
      while idx_differ < len(event_schedule):
        if event_schedule[idx_sweep_line][0] != event_schedule[idx_differ][0]:
          break
        idx_differ += 1
      */

      std::for_each(it_sweep_line, it_differ,
                    [&](const Event &e)
                    {
                      auto pos = std::get<1>(e);
                      if ((pos & 1) == 0) // even / 'sequence begin'
                      {
                        status.insert(pos >> 1);
                      } else // odd / 'sequence end'
                      {
                        status.erase((pos-1) >> 1);
                      }
                    });
      /*
      for idx_intersect in range(idx_sweep_line, idx_differ):
        if event_schedule[idx_intersect][2] == 'begin':
          sweep_line_status.add(event_schedule[idx_intersect][1])
        else:
          sweep_line_status.remove(event_schedule[idx_intersect][1])
      */

      auto first_id = *status.begin();
      for (auto id : status)
      {
        if (elements[first_id].superior_weight(elements[id]) or good[id])
        {
          break;
        }
        good.set(id);
      }
      it_sweep_line = it_differ;
    }
  }
  bool operator[](size_t idx) const noexcept
  {
    assert(idx < good.size());
    return good[idx];
  }
  [[nodiscard]] size_t count(void) const noexcept { return good.count(); }
};
#endif

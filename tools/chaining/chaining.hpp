#ifndef CHAINING_HPP
#define CHAINING_HPP
#include "utilities/constexpr_if.hpp"
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iterator>
#include <set>
#include <vector>

template <typename ElementClass, bool local = false>
class Chain
{
  private:
  static constexpr const size_t undef = ~static_cast<size_t>(0);
  std::vector<size_t> precursors;
  std::vector<uint32_t> scores;

  struct ElementValue
  {
    uint32_t value;
    explicit ElementValue(uint32_t _value)
    : value(_value)
    {}
  };

  class CmpIds
  {
    private:
    const std::vector<ElementClass> &elements;
    const std::vector<uint32_t> &scores;

    public:
    using is_transparent = std::true_type;
    CmpIds(const std::vector<ElementClass> &_elements,
           const std::vector<uint32_t> &_scores)
      : elements(_elements)
      , scores(_scores)
    {}
    bool operator ()(const size_t id1, const size_t id2) const
    {
      return elements[id1].inferior_secondary_endpos_tie_unequal(elements[id2]);
    }
    bool operator ()(const size_t id, const ElementValue &element) const
    {
      return elements[id].inferior_secondary_endpos(element.value);
    }
  };

  struct Iterator
  {
    using iterator_category = std::forward_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using value_type        = size_t;
    using pointer           = value_type*;
    using reference         = value_type&;

    Iterator(const std::vector<size_t> &_precursors,
             size_t _current)
      : precursors(_precursors)
      , current(_current)
    {}
    value_type operator*() const
    {
      return current;
    }
    Iterator& operator++() // Prefix increment
    {
      if (current != undef)
      {
        current = precursors[current];
      }
      return *this;
    }
    Iterator operator++(int) // Postfix increment
    {
      Iterator tmp = *this;
      if (current != undef)
      {
        current = precursors[current];
      }
      return tmp;
    }
    bool operator == (const Iterator& other) const
    {
      return current == other.current;
    }
    bool operator != (const Iterator& other) const
    {
      return not (*this == other);
    }

    private:
    const std::vector<size_t> &precursors;
    size_t current;
  };

  public:
  explicit Chain(const std::vector<ElementClass> &elements)
  : precursors(elements.size(), undef)
  , scores(elements.size(), 0)
  {
    // BST of element ids
    CmpIds cmp_ids(elements,scores);
    std::set<size_t, CmpIds> status(cmp_ids);

    auto prio = [&] (const size_t id)
                    { return scores[id] - elements[id].gap_score(); };

    // vector of pairs <position, id incl. event type>
    // id < 0 for 'sequence start' and > 0 for 'sequence end'
    using Event = std::pair<uint32_t,int32_t>;
    std::vector<Event> event_schedule{};
    event_schedule.reserve(2 * elements.size());

    assert(((elements.size()-1) << 1) + 1 <= UINT32_MAX);
    for (int32_t idx = 0; idx < static_cast<int32_t>(elements.size()); idx++)
    {
      event_schedule.emplace_back(elements[idx].primary_startpos_get(),
                                  -idx - 1);
      event_schedule.emplace_back(elements[idx].primary_endpos_get(),
                                  idx + 1);
    }
    std::ranges::sort(event_schedule);

    for (auto & it_sweep_line : event_schedule)
    {
      auto current = std::get<1>(it_sweep_line);

      if (current < 0) // process start point
      {
        current = -current - 1;
        uint32_t weight = elements[current].weight_get();

        auto pred_it = status.lower_bound(ElementValue(elements[current]
                                          .secondary_startpos_get()));
        if (pred_it == status.begin()) // no precursor
        {
          scores[current] = weight;
          precursors[current] = undef;
        }
        else
        {
          auto previous = *std::prev(pred_it);
          if constexpr (local)
          {
            int32_t score = static_cast<int32_t>(scores[previous] + weight)
                            - elements[current].gap_score(elements[previous]);
            if (score > static_cast<int32_t>(weight))
            {
              scores[current] = static_cast<uint32_t>(score);
              precursors[current] = previous;
            }
            else
            {
              scores[current] = weight;
              precursors[current] = undef;
            }
          }
          else
          {
            scores[current] = scores[previous] + weight;
            precursors[current] = previous;
          }
        }
      }

      else // process end point
      {
        current -= 1;

        auto pred_it = status.lower_bound(current);
        if (pred_it == status.begin() or
            constexpr_if<local>(prio(current), scores[current])
             > constexpr_if<local>(prio(*std::prev(pred_it)),
                                   scores[*std::prev(pred_it)]))
        {
          status.insert(current);

          auto succ_it = status.upper_bound(current);
          while (succ_it != status.end() and
                 constexpr_if<local>(prio(current), scores[current])
                  > constexpr_if<local>(prio(*succ_it), scores[*succ_it]))
          {
            auto tmp = succ_it;
            succ_it = status.upper_bound(*succ_it);
            status.erase(tmp);
          }
        }
      }
    }
  }
  [[nodiscard]] size_t precursor(size_t idx) const noexcept
  {
    assert(idx < precursors.size());
    return precursors[idx];
  }
  [[nodiscard]] uint32_t score(size_t idx) const noexcept
  {
    assert(idx < scores.size());
    return scores[idx];
  }
  [[nodiscard]] uint32_t score(void) const noexcept
  {
    return *std::ranges::max_element(scores);
  }
  Iterator begin(void) const noexcept
  {
    auto first_elem = std::distance(scores.begin(),
                                    std::ranges::max_element(scores));
    return Iterator(precursors, first_elem);
  }
  Iterator end(void) const noexcept
  {
    return Iterator(precursors, undef);
  }
  [[nodiscard]] size_t size(void) const noexcept
  {
    size_t size = 0;
    for (auto it = begin(); it != end(); ++it)
    {
      size++;
    }
    return size;
  }
};

#endif

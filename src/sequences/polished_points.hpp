#ifndef POLISHED_POINTS_HPP
#define POLISHED_POINTS_HPP
#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cassert>
#include <cstdbool>
#include <iostream>
#include <vector>
#include "utilities/mathsupport.hpp"
#include "sequences/align_polish.hpp"

class PolishedPoint
{
  private:
  size_t distance, row, aligned_len;
  public:
  PolishedPoint(void)
    : distance(0)
    , row(0)
    , aligned_len(0)
  {}
  PolishedPoint(size_t _distance, size_t _row, size_t _aligned_len)
    : distance(_distance)
    , row(_row)
    , aligned_len(_aligned_len)
  {
    assert(_aligned_len >= _row);
  }

  bool operator == (const PolishedPoint& other) const noexcept
  {
    return distance == other.distance
           && row == other.row
           && aligned_len == other.aligned_len;
  }

  bool operator != (const PolishedPoint& other) const noexcept
  {
    return !(*this == other);
  }

  std::string to_string(void) const noexcept
  {
    return "(" +
           std::to_string(distance) + "," +
           std::to_string(row) + "," +
           std::to_string(aligned_len) + "," +
           (aligned_len == 0 ? "undef"
                             : std::to_string(error_percentage_get(distance,
                                                                   aligned_len))
           ) + ")";
  }

  size_t distance_get(void) const noexcept
  {
    return distance;
  }

  size_t aligned_len_get(void) const noexcept
  {
    return aligned_len;
  }

  size_t row_get(void) const noexcept
  {
    return row;
  }

  size_t column_get(void) const noexcept
  {
    assert(aligned_len >= row);
    return aligned_len - row;
  }
};

class PolishedPoints
{
  struct Iterator
  {
    private:
      const PolishedPoint end_point{};
      int64_t point_idx;
      const std::vector<PolishedPoint> &best_ref;
    public:
      Iterator(const std::vector<PolishedPoint> &_best_ref,int64_t _point_idx)
        : point_idx(_point_idx)
        , best_ref(_best_ref)
     {
     }
     Iterator& operator++() /* prefix increment */
     {
       assert(point_idx >= 0);
       --point_idx;
       return *this;
     }
     const PolishedPoint &operator *(void) const noexcept
     {
       assert(point_idx >= 0);
       return point_idx == 0 ? end_point : best_ref[point_idx-1];
     }
     bool operator != (const Iterator& other) const noexcept
     {
       return point_idx != other.point_idx;
     }
  };
  private:
  bool defined;
  double smallest_error_percentage;
  size_t longest_aligned_len,
         context_distance,
         context_aligned_len;
  std::vector<PolishedPoint> best;
  public:
  PolishedPoints(size_t _context_distance,size_t _context_aligned_len)
    : defined(false)
    , smallest_error_percentage(0.0)
    , longest_aligned_len(0)
    , context_distance(_context_distance)
    , context_aligned_len(_context_aligned_len)
    , best({})
  {}
  void clear(size_t _context_distance,size_t _context_aligned_len)
  {
    defined = false;
    smallest_error_percentage = 0.0;
    longest_aligned_len = 0;
    context_distance = _context_distance;
    context_aligned_len = _context_aligned_len;
    best.clear();
  }
  void clear(void)
  {
    defined = false;
    smallest_error_percentage = 0.0;
    longest_aligned_len = 0;
    best.clear();
  }
  std::string to_string(void) const noexcept
  {
    std::string s{};
    s += '(' + std::to_string(smallest_error_percentage) + "," +
               std::to_string(longest_aligned_len);
    bool first_element = true;
    for (auto &&b : best)
    {
      if (first_element)
      {
        first_element = false;
        s += ",[";
      } else
      {
        s += ",";
      }
      s += b.to_string();
    }
    s += "]";
    return s;
  }
  void add(size_t _distance, size_t _row, size_t _aligned_len)
  {
    const double this_error_percentage
      = error_percentage_get(_distance + context_distance,
                             _aligned_len + context_aligned_len);
    if (defined)
    {
      if (longest_aligned_len < _aligned_len)
      {
        longest_aligned_len = _aligned_len;
        if (this_error_percentage < smallest_error_percentage)
        {
          /* longer length and small error percentage of current match
              dominate all previous entries and so we can delete these. */
          best.clear();
          smallest_error_percentage = this_error_percentage;
        } else
        {
          size_t n_rm = 0; /* number of elements removed */
          /* remove all entries with a larger error percentage. They
             are also not longer than the current match has maximum length. */
          for (size_t read_idx = 1; read_idx < best.size(); read_idx++)
          {
            assert(read_idx > n_rm);
            const PolishedPoint &cmp_pp = best[read_idx - 1 - n_rm];
            n_rm += (error_percentage_get(cmp_pp.distance_get()
                                            + context_distance,
                                          cmp_pp.aligned_len_get()
                                            + context_aligned_len)
                     >= this_error_percentage);
            best[read_idx - n_rm] = best[read_idx];
          }
          best.resize(best.size() - n_rm);
        }
        best.push_back(PolishedPoint(_distance,_row,_aligned_len));
      } else
      {
        if (smallest_error_percentage > this_error_percentage)
        {
          smallest_error_percentage = this_error_percentage;
          size_t n_rm = 0;
          /* remove all entries with a smaller aligned length than the current
             They also have a larger error percentage and does are dominated
             by the current entry */
          for (size_t read_idx = 1; read_idx < best.size(); read_idx++)
          {
            assert(read_idx > n_rm);
            n_rm += (best[read_idx - 1 - n_rm].aligned_len_get()
                     <= _aligned_len);
            best[read_idx - n_rm] = best[read_idx];
          }
          best.resize(best.size() - n_rm);
          best.push_back(PolishedPoint(_distance,_row,_aligned_len));
        }
        /* in the else case the current match does neither improve the
           aligned length value nor the error percentage and so we can ignore
           it */
      }
    } else
    {
      defined = true;
      longest_aligned_len = _aligned_len;
      smallest_error_percentage = this_error_percentage;
      best.push_back(PolishedPoint(_distance,_row,_aligned_len));
    }
  }
  bool operator == (const PolishedPoints& other) const noexcept
  {
    if (defined && !other.defined)
    {
      return false;
    }
    if (!defined && other.defined)
    {
      return false;
    }
    if (smallest_error_percentage != other.smallest_error_percentage)
    {
      return false;
    }
    if (longest_aligned_len != other.longest_aligned_len)
    {
      return false;
    }
    if (best.size() != other.best.size())
    {
      return false;
    }
    for (size_t idx = 0; idx < best.size(); idx++)
    {
      if (best[idx] != other.best[idx])
      {
        return false;
      }
    }
    return true;
  }
  bool operator != (const PolishedPoints& other) const noexcept
  {
    return !(*this == other);
  }
  size_t longest_aligned_len_get(void) const noexcept
  {
    return longest_aligned_len;
  }
  double smallest_error_percentage_get(void) const noexcept
  {
    return smallest_error_percentage;
  }
  size_t smallest_distance_get(void) const noexcept
  {
    return (longest_aligned_len/2) * (smallest_error_percentage/100.0);
  }
  Iterator begin(void) const noexcept
  {
    return Iterator(best,best.size());
  }
  Iterator end(void) const
  {
    return Iterator(best,-1);
  }
};

/* FrontValueClass need to implement the following functions:
   - size_t row_get(size_t ulen) const noexcept;
     returns the row value

   - uint64_t match_history_get(void) const noexcept;
     return the match history

   - size_t aligned_len_get(size_t idx, size_t d, size_t ulen,
                            size_t vlen) const noexcept
     returns the sum of the length of the aligned sequences.
     Here d is the distance value of the front to which the value belongs
     and idx is the index of this entry in a array indexed from 0 to 2 * d.

  - std::string to_string(void) const noexcept;
    return a string representation of the entry, used for debugging
*/

template<class FrontValueClass>
class TrackPolishedPoints
{
  private:
    PolishedPoints *best_polished_points;
    const AlignmentPolishing &alignment_polishing;
    size_t ulen,
           vlen,
           lag_last_d_with_pp,
           last_d_with_polished_point;
    size_t weight_frontentries(size_t d,const FrontValueClass *destfront)
    {
      size_t strong_history = 0;
      for (size_t idx = 0; idx <= 2 * d; idx++)
      {
        const FrontValueClass &front = destfront[idx];
        const uint64_t match_history = front.match_history_get();
#undef SKDEBUG
#ifdef SKDEBUG
        std::cout << "evaluatefrontentry(d=" << d << "\tdiag="
                  << (static_cast<long>(idx) - static_cast<long>(d))
                  << "\trow=" << front.to_string() << ")";
#endif
        if (alignment_polishing.is_polished(match_history))
        {
          strong_history++;
          const size_t dest_row = front.row_get(ulen);
          const size_t aligned_len = front.aligned_len_get(idx,d,ulen,vlen);
          best_polished_points->add(d,dest_row,aligned_len);
#ifdef SKDEBUG
          std::cout << "*";
          std::cout << std::endl;
          std::cout << "add polishing point for d=" << d << "\tdest_row="
                    << dest_row << "\taligned_len=" << aligned_len << std::endl;
          std::cout << "current list of polished points" << std::endl;
          std::cout << best_polished_points->to_string() << std::endl;
#endif
        }
#ifdef SKDEBUG
        else
        {
          std::cout << std::endl;
        }
#endif
      }
      return strong_history;
    }
  public:
    TrackPolishedPoints(PolishedPoints *_best_polished_points,
                        const AlignmentPolishing &_alignment_polishing,
                        size_t _ulen,
                        size_t _vlen)
      : best_polished_points(_best_polished_points)
      , alignment_polishing(_alignment_polishing)
      , ulen(_ulen)
      , vlen(_vlen)
      , lag_last_d_with_pp(_alignment_polishing.lag_last_d_with_pp_get())
      , last_d_with_polished_point(0)
     {}
  size_t evaluate(size_t d,const FrontValueClass *front)
  {
    if (d == 0)
    {
      return 1;
    }
    const size_t strong_history = weight_frontentries(d,front);
    if (strong_history > 0)
    {
      last_d_with_polished_point = d;
      return 1;
    }
    assert(last_d_with_polished_point < d);
    if (d - last_d_with_polished_point <= lag_last_d_with_pp)
    {
      return 1;
    }
    return 0;
  }
};
#endif

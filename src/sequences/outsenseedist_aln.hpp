#ifndef OUTSENSEEDIST_ALN_HPP
#define OUTSENSEEDIST_ALN_HPP
#include <cassert>
#include <cstdint>
#include <utility>
#include <vector>
#include <cstdio>
#include <cmath>
#include <iostream>
#include "sequences/backreference.hpp"
#include "sequences/front_value_trace.hpp"
#include "sequences/eoplist.hpp"

class TrackEditoperations
{
  private:
  std::vector<Backreference> trace;
  [[nodiscard]] size_t max_d_get(void) const noexcept
  {
    const size_t max_d_plus1 = std::sqrt(trace.size());
    assert(max_d_plus1 > 0 && max_d_plus1 * max_d_plus1 == trace.size());
    return max_d_plus1 - 1;
  }
  public:
    TrackEditoperations(void)
      : trace({})
     {}
  size_t evaluate([[maybe_unused]] size_t d,int32_t lo_diag, int32_t hi_diag,
                  const FrontValueTrace *front)
  {
    for (int32_t idx = lo_diag; idx <= hi_diag; idx++)
    {
      trace.push_back(front[idx].backreference_get());
    }
    return size_t(1);
  }
  void show(void) const noexcept
  {
    size_t this_offset = 0;
    for (size_t d = 0; this_offset < trace.size(); d++)
    {
      std::cout << "front(d=" << d << ")";
      for (size_t idx = this_offset; idx <= this_offset + 2 * d; idx++)
      {
        std::cout << "\t" << trace[idx].to_string();
      }
      std::cout << '\n';
      this_offset += 2 * d + 1;
    }
  }
  [[nodiscard]] Eoplist traceback_one(size_t ulen, size_t vlen) const noexcept
  {
    const size_t max_d = max_d_get();
    assert(trace.size() >= max_d + 1);
#ifndef NDEBUG
    size_t u_remain = ulen;
    size_t v_remain = vlen;
#endif
    size_t front_mid = trace.size() - 1 - max_d;
    assert(front_mid + vlen >= ulen);
    int64_t diag = static_cast<int64_t>(vlen) - static_cast<int64_t>(ulen);
    Eoplist eoplist(true);
    for (size_t current_d = max_d; /* Nothing */; current_d--)
    {
      assert(diag >= -static_cast<int64_t>(current_d) &&
             std::cmp_less_equal(diag, current_d));
      const Backreference &br = trace[front_mid + diag];
      const uint32_t match_length = br.local_matchcount_get();
      if (match_length > 0)
      {
        eoplist.match_add(static_cast<size_t>(match_length));
      }
      if (br.has_mismatch())
      {
        assert(u_remain > match_length &&
               v_remain > match_length);
#ifndef NDEBUG
        u_remain -= (match_length + 1);
        v_remain -= (match_length + 1);
#endif
        eoplist.mismatch_add();
      } else
      {
        if (br.has_deletion())
        {
          diag++;
          assert(u_remain > match_length &&
                 v_remain >= match_length);
#ifndef NDEBUG
          u_remain -= (match_length + 1);
          v_remain -= match_length;
#endif
          eoplist.deletion_add();
        } else
        {
          if (br.has_insertion())
          {
            diag--;
            assert(u_remain >= match_length &&
                   v_remain > match_length);
#ifndef NDEBUG
            u_remain -= match_length;
            v_remain -= (match_length+1);
#endif
            eoplist.insertion_add();
          } else
          {
            /* Now at d=0 and h=0 */
            assert(current_d == 0);
#ifndef NDEBUG
            u_remain = u_remain <= match_length ? 0 : (u_remain - match_length);
            v_remain = v_remain <= match_length ? 0 : (v_remain - match_length);
#endif
            break;
          }
        }
      }
      assert(front_mid >= 2 * current_d);
      front_mid -= 2 * current_d;
    }
    assert(diag == 0);
    assert(u_remain == 0);
    assert(v_remain == 0);
    eoplist.reverse_end(0);
    return eoplist;
  }
  private:
  struct StackElem
  {
    size_t u_remain,
           v_remain,
           current_d,
           front_mid;
           int64_t diag;
    StackElem(size_t _u_remain,
              size_t _v_remain,
              size_t _current_d,
              size_t _front_mid,
              int64_t _diag)
      : u_remain(_u_remain)
      , v_remain(_v_remain)
      , current_d(_current_d)
      , front_mid(_front_mid)
      , diag(_diag)
      {}
  };
  void expand(std::vector<StackElem> *stack,
              size_t u_remain, size_t v_remain, size_t current_d,
              size_t front_mid, int64_t diag) const noexcept
  {
    const Backreference &br = trace[front_mid + diag];
    const uint32_t match_length = br.local_matchcount_get();
    bool has_pushed = true;
    if (br.has_mismatch())
    {
      assert(current_d > 0);
      assert(u_remain > match_length &&
             v_remain > match_length &&
             front_mid >= 2 * current_d);
      stack->emplace_back(u_remain - (match_length + 1),
                          v_remain - (match_length + 1),
                          current_d - 1,
                          front_mid - 2 * current_d,
                          diag);
      has_pushed = true;
    }
    if (br.has_deletion())
    {
      assert(current_d > 0);
      assert(u_remain > match_length &&
             v_remain >= match_length &&
             front_mid >= 2 * current_d);
      stack->emplace_back(u_remain - (match_length + 1),
                          v_remain - match_length,
                          current_d - 1,
                          front_mid - 2 * current_d,
                          diag + 1);
      has_pushed = true;
    }
    if (br.has_insertion())
    {
      assert(current_d > 0);
      assert(u_remain >= match_length &&
             v_remain > match_length &&
             front_mid >= 2 * current_d);
      stack->emplace_back(u_remain - match_length,
                          v_remain - (match_length + 1),
                          current_d - 1,
                          front_mid - 2 * current_d,
                          diag - 1);
      has_pushed = true;
    }
    if (!has_pushed)
    {
      assert(current_d == 0);
      assert(diag == 0);
      assert(u_remain == match_length);
      assert(v_remain == match_length);
    }
  }
  public:
  void traceback_all(size_t ulen, size_t vlen) const noexcept
  {
    const size_t max_d = max_d_get();
    assert(trace.size() >= max_d + 1);
    const size_t front_mid = trace.size() - 1 - max_d;
    assert(front_mid + vlen >= ulen);
    const int64_t diag = static_cast<int64_t>(vlen) -
                         static_cast<int64_t>(ulen);
    std::vector<StackElem> stack{};
    expand(&stack,ulen,vlen,max_d,front_mid,diag);
    while (not stack.empty())
    {
      const StackElem se = stack[stack.size() - 1];
      stack.pop_back();
      expand(&stack,se.u_remain,se.v_remain,
             se.current_d,se.front_mid,se.diag);
    }
  }
};
#endif

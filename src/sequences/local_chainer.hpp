#ifndef LOCAL_CHAINER_HPP
#define LOCAL_CHAINER_HPP
#include <cassert>
#include <utility>
#include <tuple>
#include <vector>
#include <cstdint>
#include <climits>
#include <cstdlib>
#include <algorithm>
#include "utilities/unused.hpp"

template<typename GapType,typename ScoreType,typename PredecessorType>
class LocalChainElemInfo
{
  static constexpr const GapType GapTypeMax = ~static_cast<GapType>(0);
  static constexpr const PredecessorType PredecessorTypeMax
    = ~static_cast<PredecessorType>(0);
  static constexpr const PredecessorType first_bit
    = static_cast<PredecessorType>(1) <<
      (sizeof(PredecessorType) * CHAR_BIT - 1);
  using LocalChainElemInfoType = struct {GapType ref_gap_length,
                                                 query_gap_length;
                                         ScoreType score;
                                         PredecessorType predecessor;
                                         bool is_referenced;};
  LocalChainElemInfoType *chain_elem_info_array;
  size_t allocated, nextfree;
  public:
  LocalChainElemInfo(void)
    : chain_elem_info_array(nullptr)
    , allocated(0)
    , nextfree(0)
  {}
  ~LocalChainElemInfo(void)
  {
    free(chain_elem_info_array);
  }
  void free_space(void)
  {
    free(chain_elem_info_array);
    chain_elem_info_array = nullptr;
    nextfree = allocated = 0;
  }
  ScoreType score_get(size_t idx) const noexcept
  {
    assert(idx < nextfree);
    return chain_elem_info_array[idx].score;
  }
  ScoreType total_score_get(void) const noexcept
  {
    ScoreType score_sum = 0;
    for (size_t idx = 0; idx < nextfree; idx++)
    {
      score_sum += score_get(idx);
    }
    return score_sum;
  }
  bool is_referenced_get(size_t idx) const noexcept
  {
    assert(idx < nextfree);
    return chain_elem_info_array[idx].is_referenced;
  }
  void is_referenced_set(size_t idx)
  {
    assert(idx < nextfree);
    chain_elem_info_array[idx].is_referenced = true;
  }
  void is_referenced_unset(size_t idx)
  {
    assert(idx < nextfree);
    chain_elem_info_array[idx].is_referenced = false;
  }
  PredecessorType predecessor_get(size_t idx) const noexcept
  {
    assert(idx < nextfree &&
           !(chain_elem_info_array[idx].predecessor & first_bit));
    return chain_elem_info_array[idx].predecessor;
  }
  void score_set(size_t idx,ScoreType this_score)
  {
    assert(idx < nextfree);
    chain_elem_info_array[idx].score = this_score;
  }
  void chain_elem_init(size_t idx)
  {
    assert(idx < nextfree &&
           idx <= PredecessorTypeMax &&
           !(static_cast<PredecessorType>(idx) & first_bit));
    chain_elem_info_array[idx].ref_gap_length = 0;
    chain_elem_info_array[idx].query_gap_length = 0;
    chain_elem_info_array[idx].predecessor = static_cast<PredecessorType>(idx);
    chain_elem_info_array[idx].is_referenced = false;
  }
  void predecessor_set(GTTL_DEBUG_USED bool upwards_chaining,
                       size_t idx,PredecessorType p,
                       uint64_t ref_gap_length,uint64_t query_gap_length)
  {
    assert(idx < nextfree &&
           ref_gap_length <= GapTypeMax &&
           query_gap_length <= GapTypeMax);
    assert((upwards_chaining && idx > static_cast<size_t>(p)) ||
           (!upwards_chaining && idx < static_cast<size_t>(p)));
    chain_elem_info_array[idx].predecessor = p;
    chain_elem_info_array[idx].ref_gap_length
      = static_cast<GapType>(ref_gap_length);
    chain_elem_info_array[idx].query_gap_length
      = static_cast<GapType>(query_gap_length);
  }
  void marked_set(size_t idx)
  {
    assert(idx < nextfree);
    chain_elem_info_array[idx].predecessor |= first_bit;
  }
  void marked_unset(size_t idx)
  {
    assert(idx < nextfree);
    chain_elem_info_array[idx].predecessor &= ~first_bit;
  }
  bool is_marked(size_t idx) const noexcept
  {
    assert(idx < nextfree);
    return static_cast<bool>(chain_elem_info_array[idx].predecessor
                             & first_bit);
  }
  GapType ref_gap_length_get(size_t idx) const noexcept
  {
    assert(idx < nextfree);
    return chain_elem_info_array[idx].ref_gap_length;
  }
  GapType query_gap_length_get(size_t idx) const noexcept
  {
    assert(idx < nextfree);
    return chain_elem_info_array[idx].query_gap_length;
  }
  void resize(size_t segment_length)
  {
    nextfree = segment_length;
    if (allocated < segment_length)
    {
      chain_elem_info_array
        = static_cast<LocalChainElemInfoType *>
                     (realloc(chain_elem_info_array,
                              segment_length * sizeof *chain_elem_info_array));
      allocated = segment_length;
    }
  }
};

template<typename ScoreType>
class LocalChainAttributes
{
  size_t ref_endpos,
         ref_match_length,
         query_endpos,
         query_match_length;
  ScoreType score;
  size_t from_element,
         chain_len;
  public:
  LocalChainAttributes(size_t _ref_endpos,
                       size_t _ref_match_length,
                       size_t _query_endpos,
                       size_t _query_match_length,
                       ScoreType _score,
                       size_t _from_element,
                       size_t _chain_len)
    : ref_endpos(_ref_endpos)
    , ref_match_length(_ref_match_length)
    , query_endpos(_query_endpos)
    , query_match_length(_query_match_length)
    , score(_score)
    , from_element(_from_element)
    , chain_len(_chain_len)
  {}
  ScoreType score_get(void) const noexcept
  {
    return score;
  }
  bool self_overlap(void) const noexcept
  {
    return ref_endpos >= query_endpos - query_match_length + 1;
  }
  size_t from_element_get(void) const noexcept
  {
    return from_element;
  }
  size_t size(void) const noexcept
  {
    return chain_len;
  }
  size_t ref_endpos_get(void) const noexcept
  {
    return ref_endpos;
  }
  size_t ref_match_length_get(void) const noexcept
  {
    return ref_match_length;
  }
  size_t query_endpos_get(void) const noexcept
  {
    return query_endpos;
  }
  size_t query_match_length_get(void) const noexcept
  {
    return query_match_length;
  }
};

#ifndef NDEBUG
template<class SeedTable>
static void check_chain_element_order(const SeedTable &seed_table,
                                      size_t segment_start,
                                      size_t segment_length)
{
  assert(segment_length > 0);
  size_t prev_j_endpos = seed_table.order_endpos_get(segment_start);
  for (size_t j = 1; j < segment_length; j++)
  {
    const size_t j_endpos = seed_table.order_endpos_get(segment_start + j);
    assert(prev_j_endpos <= j_endpos);
    prev_j_endpos = j_endpos;
  }
}
#endif

template<class SeedTable,
         typename GapType,typename ScoreType,typename PredecessorType>
class LocalChainer
{
  using ThisLocalChainElemInfo
    = LocalChainElemInfo<GapType,ScoreType,PredecessorType>;
  using EndLengthPair = std::pair<PredecessorType,uint32_t>;
  static constexpr const ScoreType ScoreTypeMax = ~static_cast<ScoreType>(0);
  private:
  ScoreType gap_function(ScoreType ref_gap_length,ScoreType query_gap_length)
    const noexcept
  {
    return static_cast<ScoreType>(ref_gap_length + query_gap_length)/2;
  }
  const SeedTable &seed_table;
  size_t segment_start, /* offset for all access to seed_table */
         segment_length, /* number of elements of seed_table segment */
         max_previous;
  ThisLocalChainElemInfo chain_elem_info_fwd,
                         chain_elem_info_bck,
                         &chain_elem_info;
  std::vector<EndLengthPair> chain_ends_lengths;

  template<bool upwards_chaining>
  void local_chain_scores(ThisLocalChainElemInfo *var_chain_elem_info)
  {
    var_chain_elem_info->resize(segment_length);
    size_t j;
    int step;
    size_t end_segment;
    if constexpr (upwards_chaining)
    {
      j = 0;
      step = 1;
      end_segment = segment_length - 1;
    } else
    {
      j = segment_length - 1;
      step = -1;
      end_segment = 0;
    }
    const size_t match_length = seed_table.length_get(segment_start + j);
    var_chain_elem_info->score_set(j,static_cast<ScoreType>(match_length));
    var_chain_elem_info->chain_elem_init(j);
    j += step;
    if (segment_length > 1)
    {
      while (true)
      {
        size_t last_idx;
        if constexpr (upwards_chaining)
        {
          last_idx = (j <= max_previous ? 0 : (j - max_previous));
        } else
        {
          last_idx = (j + max_previous <= segment_length - 1)
                        ? (j + max_previous)
                        : (segment_length - 1);
        }
        assert(j < segment_length);
        const size_t j_match_length = seed_table.length_get(segment_start + j);
        ScoreType j_maxscore = static_cast<ScoreType>(j_match_length);
        var_chain_elem_info->chain_elem_init(j);
        assert(step < 0 || static_cast<size_t>(step) <= j);
        size_t i = j - step;
        while (true)
        {
          assert(i < segment_length);
          uint64_t ref_gap_length, query_gap_length;
          if constexpr (upwards_chaining)
          {
            std::tie(ref_gap_length,query_gap_length)
              = seed_table.colinear_match_pair(segment_start + i,
                                                segment_start + j);
          } else
          {
            std::tie(ref_gap_length,query_gap_length)
              = seed_table.colinear_match_pair(segment_start + j,
                                                segment_start + i);
          }
          assert(ref_gap_length <= static_cast<uint64_t>(ScoreTypeMax) &&
                 query_gap_length <= static_cast<uint64_t>(ScoreTypeMax) &&
                 ref_gap_length <= static_cast<uint64_t>(ScoreTypeMax)
                                   - query_gap_length);
          if (ref_gap_length > 0 || query_gap_length > 0)
          {
#ifdef SKDEBUG
            printf("ref_gap\t%zu\n",
                   static_cast<size_t>(ref_gap_length));
            printf("query_gap\t%zu\n",
                   static_cast<size_t>(query_gap_length));
#endif
#define DIFF_CHECK
#ifdef DIFF_CHECK
            const ScoreType diff = ref_gap_length >= query_gap_length
                                     ? ref_gap_length - query_gap_length
                                     : query_gap_length - ref_gap_length;
            if (diff < 100 ||
                static_cast<double>(diff)/std::max(ref_gap_length,
                                                   query_gap_length) <= 0.3)
#endif
            {
              const ScoreType gap_score = gap_function(ref_gap_length,
                                                       query_gap_length);
              const ScoreType scorefrom_i
                = var_chain_elem_info->score_get(i) +
                  static_cast<ScoreType>(j_match_length);
#ifdef SKDEBUG
              printf("gap_score\t%zu\n",static_cast<size_t>(gap_score));
              printf("score_get(%zu)\t%ld\n",segment_start + i,
                     static_cast<long>(var_chain_elem_info->score_get(i)));
              printf("j_match_length\t%zu\n",j_match_length);
              printf("maxscore\t%ld\n",static_cast<long>(j_maxscore));
#endif
              if (scorefrom_i > gap_score &&
                  j_maxscore < scorefrom_i - gap_score)
              {
                j_maxscore = scorefrom_i - gap_score;
                var_chain_elem_info->predecessor_set(upwards_chaining,/*unused*/
                                                     j,i,ref_gap_length,
                                                     query_gap_length);
              }
            }
          }
          if (i == last_idx)
          {
            break;
          }
          assert(step < 0 || static_cast<size_t>(step) <= i);
          i -= step;
        }
        var_chain_elem_info->score_set(j,j_maxscore);
        if (j == end_segment)
        {
          break;
        }
        j += step;
      }
    }
  }
  ThisLocalChainElemInfo &chain_elem_info_max_get(void)
  {
    local_chain_scores<true>(&chain_elem_info_fwd);
    local_chain_scores<false>(&chain_elem_info_bck);
    if (chain_elem_info_fwd.total_score_get()
        >= chain_elem_info_bck.total_score_get())
    {
      chain_elem_info_bck.free_space();
      return chain_elem_info_fwd;
    } else
    {
      chain_elem_info_fwd.free_space();
      return chain_elem_info_bck;
    }
  }
  public:
  LocalChainer(const SeedTable &_seed_table,
        size_t _segment_start,
        size_t _segment_length,
        size_t _max_previous)
    : seed_table(_seed_table)
    , segment_start(_segment_start)
    , segment_length(_segment_length)
    , max_previous(_max_previous)
    , chain_elem_info_fwd({})
    , chain_elem_info_bck({})
    , chain_elem_info(chain_elem_info_max_get())
    , chain_ends_lengths({})
  {
#ifndef NDEBUG
    check_chain_element_order<SeedTable>(seed_table,segment_start,
                                          segment_length);
#endif
    std::vector<PredecessorType> chain_ends{};
    const bool upwards_chaining
      = chain_elem_info_fwd.total_score_get()
        >= chain_elem_info_bck.total_score_get();
    /* As we reuse the chain_elem_info struct, we assume that
       no references are set. But we better check this. */
#ifndef NDEBUG
    for (size_t idx = 0; idx < segment_length; idx++)
    {
      assert(!chain_elem_info.is_referenced_get(idx));
    }
#endif
    /* We need to known if an element is referenced by another
       as part of a chain which of course does not end with this element.
       So whenever, the element at idx refers to a different element
       p we track this in p */
    for (size_t idx = 0; idx < segment_length; idx++)
    {
      const PredecessorType p = chain_elem_info.predecessor_get(idx);
      if (static_cast<size_t>(p) != idx)
      {
        chain_elem_info.is_referenced_set(p);
      }
    }
    /* We next collect all elements which are not referenced, so that
       they could be ends of chains. */
    for (size_t idx = 0; idx < segment_length; idx++)
    {
      if (chain_elem_info.is_referenced_get(idx))
      {
        chain_elem_info.is_referenced_unset(idx);
      } else
      {
        chain_ends.push_back(static_cast<PredecessorType>(idx));
      }
    }
    /* Next we sort the ends of the chains in descending order of their
       score. For two elements with the same score are sorted in
       ascending order of the index. */
    std::sort(chain_ends.begin(),chain_ends.end(),
              [&](PredecessorType a, PredecessorType b)
              {
                return (chain_elem_info.score_get(static_cast<size_t>(a)) >
                        chain_elem_info.score_get(static_cast<size_t>(b))) ||
                       (chain_elem_info.score_get(static_cast<size_t>(a)) ==
                        chain_elem_info.score_get(static_cast<size_t>(b)) &&
                        a < b);
              });
    /* We now process the chain end elements by descending order of their
       score. For each such element with trace back the
       chaia, mark each element until the element refers to itself (so
       that the chain cannot be extended further or a marked element
       is detected. We track the score, overwrite the previously stored
       score (as we are only interested in scores of non-overlapping chaings)
       and store the chain length with the index. */
    for (auto ce : chain_ends)
    {
      size_t j = static_cast<size_t>(ce);
      int64_t score_sum
        = static_cast<int64_t>(seed_table.length_get(segment_start + j));
      uint32_t chain_len = 1;
      while (true)
      {
        const PredecessorType i = chain_elem_info.predecessor_get(j);
        chain_elem_info.marked_set(j);
        if (i == j || chain_elem_info.is_marked(i))
        {
          break;
        }
        const ScoreType gap_score
          = gap_function(chain_elem_info.ref_gap_length_get(j),
                         chain_elem_info.query_gap_length_get(j));
        score_sum -= static_cast<int64_t>(gap_score);
        score_sum
          += static_cast<int64_t>(seed_table.length_get(segment_start + i));
        chain_len++;
        j = i;
      }
      assert(score_sum <= static_cast<int64_t>(chain_elem_info.score_get(ce)));
      chain_elem_info.score_set(ce,static_cast<ScoreType>(std::max(int64_t(0),
                                                                   score_sum)));
      chain_ends_lengths.push_back(std::make_pair(ce,chain_len));
    }
    /* We now do not need the markings anymore and so we unset it */
    for (size_t idx = 0; idx < segment_length; idx++)
    {
      chain_elem_info.marked_unset(idx);
    }
    if (!upwards_chaining)
    {
      /* In this case we have constructed the chains by processing the
         matches in reverse order, so that the predecessor have a larger
         index that the origin. We want to further process the chains
         in the same way as before, so we just exchange origin and
         predecessor. Before we do this, we cheeck that each element
         belongs to exactly one non-overlapping chain. */
#ifndef NDEBUG
      std::vector<size_t> belongs_to(segment_length,segment_length);
      for (size_t idx = 0; idx < chain_ends_lengths.size(); idx++)
      {
        auto ce = chain_ends_lengths[idx];
        auto chain_len = std::get<1>(ce);
        size_t j = std::get<0>(ce);
        belongs_to[j] = j;

        if (chain_len == 1)
        {
          continue;
        }
        PredecessorType i = chain_elem_info.predecessor_get(j);
        for (size_t chain_elem_idx = 1; chain_elem_idx < chain_len && i != j;
             chain_elem_idx++)
        {
          assert(belongs_to[i] == segment_length);
          belongs_to[i] = std::get<0>(ce);
          j = i;
          i = chain_elem_info.predecessor_get(j);
        }
      }
#endif
      for (size_t idx = 0; idx < chain_ends_lengths.size(); idx++)
      {
        auto ce = chain_ends_lengths[idx];
        auto chain_len = std::get<1>(ce);

        if (chain_len == 1)
        {
          continue;
        }
        size_t j = std::get<0>(ce);
        PredecessorType i = chain_elem_info.predecessor_get(j);
        assert(i >= j);
        auto j_ref_gap_length = chain_elem_info.ref_gap_length_get(j);
        auto j_query_gap_length = chain_elem_info.query_gap_length_get(j);
        assert (belongs_to[i] == j);
        for (size_t chain_elem_idx = 1; chain_elem_idx < chain_len && i != j;
             chain_elem_idx++)
        {
          assert (belongs_to[i] == std::get<0>(ce));
          const PredecessorType next_pred = chain_elem_info.predecessor_get(i);
          auto i_ref_gap_length = chain_elem_info.ref_gap_length_get(i);
          auto i_query_gap_length = chain_elem_info.query_gap_length_get(i);
          /* We here consider the case of downwards_chaining, i.e.
             upwards_chaining = false and we want to convert it such that
             the condition of upwards_chaining are true, so we
             supply !upwards_training = true as first argument */
          chain_elem_info.predecessor_set(!upwards_chaining,i,j,
                                          j_ref_gap_length,
                                          j_query_gap_length);
          j = i;
          j_ref_gap_length = i_ref_gap_length;
          j_query_gap_length = i_query_gap_length;
          i = next_pred;
        }
        chain_ends_lengths[idx] = std::make_pair<PredecessorType,uint32_t>
                                             (static_cast<PredecessorType>(i),
                                              static_cast<uint32_t>(chain_len));
        chain_elem_info.score_set(i,
                                  chain_elem_info.score_get(std::get<0>(ce)));
        chain_elem_info.chain_elem_init(std::get<0>(ce));
      }
    }
    /* The following now sorts the ends of the chains (together with the
       length of the chains) again by their previously updated score. */
    std::sort(chain_ends_lengths.begin(),chain_ends_lengths.end(),
              [&](EndLengthPair a, EndLengthPair b)
              {
                const PredecessorType pa = std::get<0>(a);
                const PredecessorType pb = std::get<0>(b);

                return (chain_elem_info.score_get(static_cast<size_t>(pa)) >
                        chain_elem_info.score_get(static_cast<size_t>(pb))) ||
                       (chain_elem_info.score_get(static_cast<size_t>(pa)) ==
                        chain_elem_info.score_get(static_cast<size_t>(pb)) &&
                        pa < pb);

              });
  }
  /* The following is only used when chain_closed or chain_elements display
     option is set. */
  auto local_chain_get(size_t from_element,uint32_t chain_len) const noexcept
  {
    std::vector<PredecessorType> chain_elements{};
    PredecessorType j = static_cast<PredecessorType>(from_element);
    assert(chain_len >= 1);
    chain_elements.push_back(j);
    while (chain_elements.size() < chain_len)
    {
      j = chain_elem_info.predecessor_get(static_cast<size_t>(j));
      chain_elements.push_back(j);
    }
    return chain_elements;
  }
  GapType ref_gap_length_get(size_t idx) const noexcept
  {
    return chain_elem_info.ref_gap_length_get(idx);
  }
  GapType query_gap_length_get(size_t idx) const noexcept
  {
    return chain_elem_info.query_gap_length_get(idx);
  }
  PredecessorType predecessor_get(size_t idx) const noexcept
  {
    return chain_elem_info.predecessor_get(idx);
  }
  class Iterator
  {
    const SeedTable &seed_table;
    size_t segment_start;
    const ThisLocalChainElemInfo &chain_elem_info_ref;
    const std::vector<EndLengthPair> &chain_ends_lengths_ref;
    size_t current_idx;
#ifndef NDEBUG
    ScoreType previous_score;
#endif

    size_t first_in_chain_idx_get(size_t from_element,
                                  size_t chain_len) const noexcept
    {
      assert(chain_len >= 2);
      size_t j = from_element;
      for (size_t current_length = 1; current_length < chain_len;
           current_length++)
      {
        j = static_cast<size_t>(chain_elem_info_ref.predecessor_get(j));
      }
      return j;
    }
    public:
    Iterator(const SeedTable &_seed_table,
             const ThisLocalChainElemInfo &_chain_elem_info_ref,
             const std::vector<EndLengthPair> &_chain_ends_lengths_ref,
             size_t _segment_start,
             size_t _current_idx)
     : seed_table(_seed_table)
     , segment_start(_segment_start)
     , chain_elem_info_ref(_chain_elem_info_ref)
     , chain_ends_lengths_ref(_chain_ends_lengths_ref)
     , current_idx(_current_idx)
#ifndef NDEBUG
     , previous_score(~static_cast<ScoreType>(0))
#endif
    {}
    LocalChainAttributes<ScoreType> operator*(void)
#ifdef NDEBUG
      const noexcept
#endif
    {
      assert(current_idx < chain_ends_lengths_ref.size());
      const EndLengthPair chain_end_length
        = chain_ends_lengths_ref[current_idx];
      auto chain_end = std::get<0>(chain_end_length);
#ifndef NDEBUG
      assert(chain_elem_info_ref.score_get(chain_end) <= previous_score);
      previous_score = chain_elem_info_ref.score_get(chain_end);
#endif
      auto chain_len = std::get<1>(chain_end_length);
      if (chain_len == 1)
      {
        const size_t ce = segment_start + chain_end;
        return LocalChainAttributes<ScoreType>
                          (seed_table.ref_endpos_get(ce),
                           seed_table.length_get(ce),
                           seed_table.query_endpos_get(ce),
                           seed_table.length_get(ce),
                           chain_elem_info_ref.score_get(chain_end),
                           chain_end,
                           chain_len);
      }
      const size_t first_chain_idx
        = first_in_chain_idx_get(chain_end,chain_len);
      const uint64_t first_match_length
        = seed_table.length_get(segment_start + first_chain_idx);
      const uint64_t ref_chain_start
        = seed_table.ref_endpos_get(segment_start + first_chain_idx)
          + 1 - first_match_length;
      const uint64_t query_chain_start
        = seed_table.query_endpos_get(segment_start + first_chain_idx)
          + 1 - first_match_length;
      const uint64_t ref_chain_end
        = seed_table.ref_endpos_get(segment_start + chain_end);
      const uint64_t query_chain_end
        = seed_table.query_endpos_get(segment_start + chain_end);
      return LocalChainAttributes<ScoreType>
                                 (static_cast<size_t>(ref_chain_end),
                                  ref_chain_end - ref_chain_start + 1,
                                  static_cast<size_t>(query_chain_end),
                                  query_chain_end - query_chain_start + 1,
                                  chain_elem_info_ref.score_get(chain_end),
                                  chain_end,
                                  chain_len);
    }
    Iterator &operator++() /* prefix increment*/
    {
      current_idx++;
      return *this;
    }
    bool operator != (const Iterator& other) const noexcept
    {
      return current_idx != other.current_idx;
    }
  };
  Iterator begin(void) const noexcept
  {
    return Iterator(seed_table,chain_elem_info,chain_ends_lengths,
                    segment_start,0);
  }
  Iterator end(void) const noexcept
  {
    return Iterator(seed_table,chain_elem_info,chain_ends_lengths,
                    segment_start,chain_ends_lengths.size());
  }
};
#endif

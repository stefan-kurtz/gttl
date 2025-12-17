#ifndef EOPLIST_HPP
#define EOPLIST_HPP
#include <cstddef>
#include <cstdint>
#include <cctype>
#include <cassert>
#include <stdexcept>
#include <utility>
#include <vector>
#include <string>
#include <algorithm>

enum EopType : uint8_t
{
  DeletionOp = 0,
  InsertionOp,
  MismatchOp,
  MatchOp,
  UndefinedOp
};

class CigarOperator
{
 public:
  static constexpr const char deletion_char = 'D';
  static constexpr const char insertion_char = 'I';
  static constexpr const char match_char = '=';
  static constexpr const char mismatch_char = 'X';
  static constexpr const char replacement_char = 'M';
  private:
  static constexpr const char eoplist_pretty_print_dm[] = {deletion_char,
                                                           insertion_char,
                                                           mismatch_char,
                                                           match_char};
  static constexpr const char eoplist_pretty_print[] = {deletion_char,
                                                        insertion_char,
                                                        replacement_char,
                                                        replacement_char};

  public:
  EopType edit_operation;
  size_t iteration;
  CigarOperator(void)
    : edit_operation(UndefinedOp)
    , iteration(0)
  {}
  [[nodiscard]] char to_char(bool distinguish_mismatch_match) const noexcept
  {
    return distinguish_mismatch_match
             ? eoplist_pretty_print_dm[static_cast<int>(edit_operation)]
             : eoplist_pretty_print[static_cast<int>(edit_operation)];
  }
  [[nodiscard]] std::string
  to_string(bool distinguish_mismatch_match) const noexcept
  {
    std::string s{};
    assert(edit_operation < UndefinedOp);
    s += std::to_string(iteration);
    s += this->to_char(distinguish_mismatch_match);
    return s;
  }
};

class Eoplist
{
  static constexpr const uint8_t ft_eopcode_maxmatches = uint8_t(253);
  static constexpr const uint8_t ft_eopcode_mismatch = uint8_t(253);
  static constexpr const uint8_t ft_eopcode_deletion = uint8_t(254);
  static constexpr const uint8_t ft_eopcode_insertion = uint8_t(255);

  [[nodiscard]] bool eopcode_is_match(uint8_t eopcode) const noexcept
  {
    return eopcode < ft_eopcode_mismatch;
  }
  [[nodiscard]] bool eopcode_is_mismatch(uint8_t eopcode) const noexcept
  {
    return eopcode == ft_eopcode_mismatch;
  }
  [[nodiscard]] bool eopcode_is_deletion(uint8_t eopcode) const noexcept
  {
    return eopcode == ft_eopcode_deletion;
  }
  [[nodiscard]] bool eopcode_is_insertion(uint8_t eopcode) const noexcept
  {
    return eopcode == ft_eopcode_insertion;
  }
  struct Iterator
  {
    private:
    const std::vector<uint8_t> &eoplist_ref;
    bool distinguish_mismatch_match,
         exhausted;
    size_t eop_idx;
    CigarOperator cigar_operator;
    void cigar_operator_next(void)
    {
      cigar_operator.edit_operation = UndefinedOp;
      cigar_operator.iteration = 0;
      bool stop = false;
      while (!stop && eop_idx < eoplist_ref.size())
      {
        const uint8_t this_eop = eoplist_ref[eop_idx];
        if (cigar_operator.iteration > 0)
        {
          assert(cigar_operator.edit_operation != UndefinedOp);
          switch(this_eop)
          {
#define EXTEND_CIGAR_OP_OF_SAME_TYPE(EOPTYPE)\
            if (cigar_operator.edit_operation == EOPTYPE)\
            {\
              cigar_operator.iteration++;\
              eop_idx++;\
            } else\
            {\
              stop = true;\
            }
            case ft_eopcode_deletion:
              EXTEND_CIGAR_OP_OF_SAME_TYPE(DeletionOp);
              break;
            case ft_eopcode_insertion:
              EXTEND_CIGAR_OP_OF_SAME_TYPE(InsertionOp);
              break;
            case ft_eopcode_mismatch:
              if (distinguish_mismatch_match)
              {
                EXTEND_CIGAR_OP_OF_SAME_TYPE(MismatchOp);
              } else
              {
                EXTEND_CIGAR_OP_OF_SAME_TYPE(MatchOp);
              }
              break;
            default:
              if (cigar_operator.edit_operation == MatchOp)
              {
                assert(this_eop < ft_eopcode_maxmatches);
                cigar_operator.iteration += (1 + this_eop);
                eop_idx++;
              } else
              {
                stop = true;
              }
          }
        } else
        {
          switch (this_eop)
          {
            case ft_eopcode_deletion:
              cigar_operator.edit_operation = DeletionOp;
              cigar_operator.iteration = size_t(1);
              break;
            case ft_eopcode_insertion:
              cigar_operator.edit_operation = InsertionOp;
              cigar_operator.iteration = size_t(1);
              break;
            case ft_eopcode_mismatch:
              cigar_operator.edit_operation
                = distinguish_mismatch_match ? MismatchOp : MatchOp;
              cigar_operator.iteration = size_t(1);
              break;
            default:
              cigar_operator.edit_operation = MatchOp;
              cigar_operator.iteration = 1 + this_eop;
              break;
          }
          eop_idx++;
        }
      }
    }
    public:
    Iterator(const std::vector<uint8_t> &_eoplist,
             bool _distinguish_mismatch_match,bool _exhausted)
      : eoplist_ref(_eoplist)
      , distinguish_mismatch_match(_distinguish_mismatch_match)
      , exhausted(_exhausted)
      , eop_idx(0)
      , cigar_operator({})
    {
      assert(eop_idx < eoplist_ref.size());
      if (!exhausted)
      {
        cigar_operator_next();
      }
    }
    const CigarOperator &operator*(void) const noexcept
    {
      return cigar_operator;
    }
    bool operator != (const Iterator& other) const noexcept
    {
      return exhausted != other.exhausted;
    }
    Iterator& operator++() /* prefix increment*/
    {
      assert(!exhausted);
      if (eop_idx >= eoplist_ref.size())
      {
        exhausted = true;
      } else
      {
        cigar_operator_next();
      }
      return *this;
    }
  };
  std::vector<uint8_t> eoplist;
  size_t counter_for_matches,
         counter_for_mismatches,
         counter_for_deletions,
         counter_for_insertions,
         counter_for_gap_opens;
  bool previous_was_gap,
       distinguish_mismatch_match;
  void indel_add(int code)
  {
    eoplist.push_back(static_cast<uint8_t>(code));
    if (!previous_was_gap)
    {
      counter_for_gap_opens++;
      previous_was_gap = true;
    }
  }
  public:
  Eoplist(bool _distinguish_mismatch_match)
    : eoplist({})
    , counter_for_matches(0)
    , counter_for_mismatches(0)
    , counter_for_deletions(0)
    , counter_for_insertions(0)
    , counter_for_gap_opens(0)
    , previous_was_gap(false)
    , distinguish_mismatch_match(_distinguish_mismatch_match)
  {}
  void reset(void)
  {
    eoplist.clear();
    counter_for_matches = 0;
    counter_for_mismatches = 0;
    counter_for_deletions = 0;
    counter_for_insertions = 0;
    counter_for_gap_opens = 0;
    previous_was_gap = false;
  }
  void match_add(size_t length)
  {
    assert(length > 0);
    counter_for_matches += length;
    while (length > 0)
    {
      if (not eoplist.empty() and
          eoplist[eoplist.size()-1] < ft_eopcode_maxmatches - 1)
      {
        if (static_cast<size_t>(eoplist[eoplist.size()-1]) + length
              < static_cast<size_t>(ft_eopcode_maxmatches))
        {
          eoplist[eoplist.size()-1] += static_cast<uint8_t>(length);
          length = 0;
        } else
        {
          assert(std::cmp_less(ft_eopcode_maxmatches -
                                 eoplist[eoplist.size()-1],
                               length));
          length = ft_eopcode_maxmatches - eoplist[eoplist.size()-1];
          /* R max */
          eoplist[eoplist.size()-1] = ft_eopcode_maxmatches - 1;
        }
      } else
      {
        if (length <= static_cast<size_t>(ft_eopcode_maxmatches))
        {
          eoplist.push_back(static_cast<uint8_t>(length - 1)); /* R length */
          length = 0;
        } else
        {
          /* R max */
          eoplist.push_back(ft_eopcode_maxmatches-1);
          length -= static_cast<size_t>(ft_eopcode_maxmatches);
        }
      }
    }
    previous_was_gap = false;
  }
  void mismatch_add(void)
  {
    eoplist.push_back(static_cast<uint8_t>(ft_eopcode_mismatch)); /* R 1 */
    counter_for_mismatches++;
    previous_was_gap = false;
  }
  void deletion_add(void)
  {
    indel_add(ft_eopcode_deletion);
    counter_for_deletions++;
  }
  void insertion_add(void)
  {
    indel_add(ft_eopcode_insertion);
    counter_for_insertions++;
  }
  [[nodiscard]] size_t size(void) const noexcept { return eoplist.size(); }
  void reverse_end(size_t firstindex)
  {
    if (firstindex + 1 >= eoplist.size())
    {
      return;
    }
    std::reverse(eoplist.begin() + firstindex,eoplist.end());
  }
  [[nodiscard]] size_t count_matches_get(void) const noexcept
  {
    return counter_for_matches;
  }
  [[nodiscard]] size_t count_mismatches_get(void) const noexcept
  {
    return counter_for_mismatches;
  }
  [[nodiscard]] size_t count_deletions_get(void) const noexcept
  {
    return counter_for_deletions;
  }
  [[nodiscard]] size_t count_insertions_get(void) const noexcept
  {
    return counter_for_insertions;
  }
  [[nodiscard]] size_t count_gap_opens_get(void) const noexcept
  {
    return counter_for_gap_opens;
  }
  [[nodiscard]] size_t aligned_len_get(void) const noexcept
  {
    return count_deletions_get() + count_insertions_get() +
           2 * (count_mismatches_get() + count_matches_get());
  }
  [[nodiscard]] size_t aligned_len_u_get(void) const noexcept
  {
    return count_deletions_get() + count_mismatches_get() + count_matches_get();
  }
  [[nodiscard]] size_t aligned_len_v_get(void) const noexcept
  {
    return count_insertions_get() + count_mismatches_get()
                                  + count_matches_get();
  }
  [[nodiscard]] size_t errors_get(void) const noexcept
  {
    return count_deletions_get() + count_insertions_get() +
           count_mismatches_get();
  }
  [[nodiscard]] double error_percentage_get(void) const noexcept
  {
    return 200.0 * static_cast<double>(errors_get())/aligned_len_get();
  }
  [[nodiscard]] Iterator begin(void) const
  {
    return Iterator(eoplist,distinguish_mismatch_match,false);
  }
  [[nodiscard]] Iterator end(void) const
  {
    return Iterator(eoplist,distinguish_mismatch_match,true);
  }
  [[nodiscard]] std::string to_string(void) const noexcept
  {
    std::string s{};
    for (auto &&co : *this)
    {
      s += co.to_string(distinguish_mismatch_match);
    }
    return s;
  }
  bool operator ==(const Eoplist &other) const noexcept
  {
    return std::equal(eoplist.begin(),eoplist.end(),other.eoplist.begin());
  }
  bool operator !=(const Eoplist &other) const noexcept
  {
    return !(*this == other);
  }
  Eoplist(bool _distinguish_mismatch_match,const std::string &cigar_string)
    : eoplist({})
    , counter_for_matches(0)
    , counter_for_mismatches(0)
    , counter_for_deletions(0)
    , counter_for_insertions(0)
    , counter_for_gap_opens(0)
    , previous_was_gap(false)
    , distinguish_mismatch_match(_distinguish_mismatch_match)
  {
    size_t iteration = 0;

    for (auto &&cc : cigar_string)
    {
      if (std::isdigit(cc))
      {
        iteration = iteration * size_t(10) + static_cast<size_t>(cc - '0');
      } else
      {
        if (cc == CigarOperator::deletion_char)
        {
          for (size_t iter = 0; iter < iteration; iter++)
          {
            deletion_add();
          }
        } else
        {
          if (cc == CigarOperator::insertion_char)
          {
            for (size_t iter = 0; iter < iteration; iter++)
            {
              insertion_add();
            }
          } else
          {
            if (cc == CigarOperator::mismatch_char)
            {
              for (size_t iter = 0; iter < iteration; iter++)
              {
                mismatch_add();
              }
            } else
            {
              if (cc == CigarOperator::match_char ||
                  cc == CigarOperator::replacement_char)
              {
                match_add(iteration);
              } else
              {
                throw std::runtime_error(
                  "illegal symbol '" + std::to_string(cc) +
                  "' in cigar string");
              }
            }
          }
        }
        iteration = 0;
      }
    }
  }
  /* A constructor to create an Eoplist from a cigar string encoding as
     used in SAM */
  template<typename CharType>
  Eoplist(bool _distinguish_mismatch_match,
          const CharType *useq, const CharType *vseq,
          const std::vector<uint32_t> &cigar_string_encoding)
    : eoplist({})
    , counter_for_matches(0)
    , counter_for_mismatches(0)
    , counter_for_deletions(0)
    , counter_for_insertions(0)
    , counter_for_gap_opens(0)
    , previous_was_gap(false)
    , distinguish_mismatch_match(_distinguish_mismatch_match)
  {
    size_t upos = 0;
    size_t vpos = 0;

    for (auto cigar_operation : cigar_string_encoding)
    {
      const uint32_t iteration             = cigar_operation >> 4;
      const uint32_t edit_operation_number = cigar_operation & 0xf;
      assert(edit_operation_number < 3);
      if (edit_operation_number == 0) /* replacement */
      {
        for (uint32_t iter = 0; iter < iteration; iter++)
        {
          if (useq[upos] == vseq[vpos])
          {
            match_add(1);
          } else
          {
            mismatch_add();
          }
          upos++;
          vpos++;
        }
      } else
      {
        if (edit_operation_number == 1)
        {
          for (uint32_t iter = 0; iter < iteration; iter++)
          {
            insertion_add();
          }
          vpos += iteration;
        } else
        {
          for (uint32_t iter = 0; iter < iteration; iter++)
          {
            deletion_add();
          }
          upos += iteration;
        }
      }
    }
  }

  template<typename ScoreType,class SeqClass_u,class SeqClass_v>
  ScoreType evaluate_score(const SeqClass_u &useq,
                           const SeqClass_u &vseq,
                           int8_t gap_open_penalty,
                           int8_t gap_extension_penalty,
                           const int8_t * const *scorematrix2D) const
  {
    size_t idx_u = 0;
    size_t idx_v = 0;
    ScoreType sum_score = 0;
#ifndef NDEBUG
    bool last_was_indel = false;
#endif

    for (const CigarOperator &co : *this)
    {
      assert(co.iteration > 0);
      if (co.edit_operation == MatchOp || co.edit_operation == MismatchOp)
      {
        for (size_t j = 0; j < co.iteration; j++)
        {
          assert(idx_u < useq.size());
          assert(idx_v < vseq.size());
          const uint8_t cc_a = useq[idx_u];
          const uint8_t cc_b = vseq[idx_v];
          sum_score += scorematrix2D[cc_a][cc_b];
          idx_u++;
          idx_v++;
        }
#ifndef NDEBUG
        last_was_indel = false;
#endif
      } else
      {
        assert ((co.edit_operation == DeletionOp ||
                 co.edit_operation == InsertionOp) and not last_was_indel);
        sum_score -= (static_cast<ScoreType>(gap_open_penalty) +
                      static_cast<ScoreType>(co.iteration) *
                      static_cast<ScoreType>(gap_extension_penalty));
#ifndef NDEBUG
        last_was_indel = true;
#endif
        if (co.edit_operation == DeletionOp)
        {
          idx_u += co.iteration;
        } else
        {
          idx_v += co.iteration;
        }
      }
    }
    assert(idx_u == useq.size() && idx_v == vseq.size());
    return sum_score;
  }

  [[nodiscard]] std::string
  cigar_string_get(bool distinguish_mismatch_match) const noexcept
  {
    std::string cigar_string{};
    for (auto &&cigar_operator : *this)
    {
      cigar_string += cigar_operator.to_string(distinguish_mismatch_match);
    }
    return cigar_string;
  }
  bool cut_off_unpolished_tail(void)
  {
    assert(not eoplist.empty());
    size_t idx = eoplist.size()-1;
    uint8_t eopcode;

    while (not eopcode_is_match(eopcode = eoplist[idx]))
    {
      if (eopcode_is_mismatch(eopcode))
      {
        assert(counter_for_mismatches > 0);
        counter_for_mismatches--;
      } else
      {
        if (eopcode_is_deletion(eopcode))
        {
          assert(counter_for_deletions > 0);
          counter_for_deletions--;
        } else
        {
          assert(eopcode_is_insertion(eopcode) && counter_for_insertions > 0);
          counter_for_insertions--;
        }
      }
      if (idx > 0)
      {
        idx--;
      } else
      {
        assert(false); /* This case implies that all eops are
                          mismatches or indels, which is not supposed to occur
                        */
        break;
      }
    }
    const size_t diff_dist = eoplist.size() - 1 - idx;
    for (size_t j = 0; j < diff_dist; j++)
    {
      eoplist.pop_back();
    }
    return diff_dist > 0;
  }
};
#endif

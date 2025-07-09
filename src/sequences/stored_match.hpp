#ifndef STORED_MATCH_HPP
#define STORED_MATCH_HPP
#include <cstddef>
#include <string>
#include <cstdint>
#include <cassert>

/* This class represents coordinates of matches for the same segment,
   i.e. for the same pair of reference and query sequence. So the
   ID of the sequences need not be stored. The coordinates of the match
   on the primary sequence, given by a start position and a length
   is used to resolve ties when comparing matches for which the
   weights are equal. The weight is a numerical value used to
   order stored matches. The order may be ascending, i.e. the larger the
   weight, the better. This order is used when the template variable
   argest_weight_best is true. The order may be descending, i.e. the
   smaller the better. This order is used when the template variable
   argest_weight_best is false.
 */

class GttlStoredMatch
{
  private:
  static constexpr const size_t max_value = UINT32_MAX;
  const uint32_t primary_startpos,
                 primary_len,
                 secondary_startpos,
                 secondary_len,
                 distance;
  const double weight;

  [[nodiscard]] double weight_get(void) const noexcept { return weight; }

  public:
  GttlStoredMatch(size_t _primary_startpos, size_t _primary_len,
                  size_t _secondary_startpos, size_t _secondary_len,
                  size_t _distance,
                  double _weight)
    : primary_startpos(_primary_startpos)
    , primary_len(_primary_len)
    , secondary_startpos(_secondary_startpos)
    , secondary_len(_secondary_len)
    , distance(_distance)
    , weight(_weight)
  {
    assert(_primary_startpos <= max_value &&
           _primary_len <= max_value &&
           _secondary_startpos <= max_value &&
           _secondary_len <= max_value &&
           _distance <= max_value);
  }
  [[nodiscard]] uint32_t primary_startpos_get(void) const noexcept
  {
    return primary_startpos;
  }
  [[nodiscard]] uint32_t primary_len_get(void) const noexcept
  {
    return primary_len;
  }
  [[nodiscard]] uint32_t secondary_startpos_get(void) const noexcept
  {
    return secondary_startpos;
  }
  [[nodiscard]] uint32_t secondary_len_get(void) const noexcept
  {
    return secondary_len;
  }
  [[nodiscard]] uint32_t distance_get(void) const noexcept { return distance; }
  /* the == and != operators are needed to check the consistency of the filter
     results in verify_filer */
  bool operator==(const GttlStoredMatch& other) const noexcept
  {
    return primary_startpos == other.primary_startpos and
           primary_len == other.primary_len and
           secondary_startpos == other.secondary_startpos and
           secondary_len == other.secondary_len and
           distance == other.distance;
  }
  bool operator != (const GttlStoredMatch& other) const noexcept
  {
    return not (*this == other);
  }

  /* compare by weight only */
  [[nodiscard]] bool
  superior_weight(const GttlStoredMatch &other) const noexcept
  {
    return this->weight_get() > other.weight_get();
  }

  /* compare weights and if a tie, compare by primary startpos */
  [[nodiscard]] bool superior_weight_tie_primary_startpos(
                               const GttlStoredMatch &other) const noexcept
  {
      return this->superior_weight(other) or
             (this->weight_get() == other.weight_get() and
              this->primary_startpos_get() > other.primary_startpos_get());
  }

  [[nodiscard]] std::string to_string(void) const noexcept
  {
    return std::to_string(primary_startpos) + "," +
           std::to_string(primary_len) + "," +
           std::to_string(secondary_startpos) + "," +
           std::to_string(secondary_len) + "," +
           std::to_string(weight);
  }
};
#endif

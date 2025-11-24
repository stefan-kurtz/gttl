#ifndef ORDERED_INTEGER_SEQUENCE_HPP
#define ORDERED_INTEGER_SEQUENCE_HPP
#include <climits>
#include <cstdint>
#include <iostream>
#include <cassert>
#include <cstddef>
#include <cstdlib>

template <typename Basetype>
class OrderedIntegerSequence
{
  static constexpr const size_t undefined_element = ~size_t(0);
  Basetype *elements;
  int logsectionsize;
  size_t nextfree,
         maxelement,
         currentsectionnum,
         numofsections,
         nofelements,
         previouselem,
         *sectionstart;

#ifndef NDEBUG
  [[nodiscard]]
  size_t linearsearch_sec_idx_largest_leq(size_t idx) const noexcept
  {
    size_t section = 0;
    for (section = 0; section < numofsections && idx >= sectionstart[section];
         section++) /* Nothing */ ;
    assert(section > 0);
    return section - 1;
  }
#endif
  public:
  [[nodiscard]]
  size_t backward_sec_idx_largest_leq(size_t idx,size_t section_upperbound)
    const noexcept
  {
    const size_t *ssptr = sectionstart + section_upperbound;
    assert(section_upperbound < numofsections);
    while (idx < *ssptr)
    {
      assert(ssptr > sectionstart);
      ssptr--;
    }
    assert(ssptr >= sectionstart && idx >= *ssptr);
    return static_cast<size_t>(ssptr - sectionstart);
  }

  private:
  [[nodiscard]]
  size_t binarysearch_sec_idx_largest_leq(size_t idx) const noexcept
  {
    const size_t *secstart_begin = sectionstart;
    const size_t *secstart_end = sectionstart + numofsections - 1;
    const size_t *found = nullptr;
    while (secstart_begin <= secstart_end)
    {
      const size_t *const midptr = secstart_begin
                                   + static_cast<size_t>(secstart_end
                                                         - secstart_begin)
                                                                / 2;
      if (idx < *midptr)
      {
        secstart_end = midptr - 1;
      } else
      {
        found = midptr;
        secstart_begin = midptr + 1;
      }
    }
    assert(found != nullptr);
    while (found + 1 < sectionstart + numofsections && *(found+1) <= idx)
    {
      found++;
    }
    return static_cast<size_t>(found - sectionstart);
  }
  bool binarysearch_is_member(const Basetype *leftptr, const Basetype *rightptr,
                              Basetype elem) const noexcept
  {

    while (leftptr <= rightptr)
    {
      const Basetype *midptr
        = leftptr + static_cast<size_t>(rightptr - leftptr)/2;
      if (elem < *midptr)
      {
        rightptr = midptr - 1;
      } else
      {
        if (elem > *midptr)
        {
          leftptr = midptr + 1;
        } else
        {
          return true;
        }
      }
    }
    return false;
  }

  size_t binarysearch_pos2seqnum(const Basetype *leftptr,
                                 const Basetype *rightptr,
                                 Basetype pos) const noexcept
  {
    assert(leftptr <= rightptr);
    if (pos < *leftptr)
    {
      return 0;
    }
    if (pos > *rightptr)
    {
      return size_t(1) + static_cast<size_t>(rightptr - leftptr);
    }
    if (pos == *leftptr)
    {
      return 0;
    }
    if (pos == *rightptr)
    {
      return static_cast<size_t>(rightptr - leftptr);
    }
    assert(pos > *leftptr && pos < *rightptr);
    const Basetype *found = nullptr;
    const Basetype *const leftorig = leftptr;
    while (leftptr <= rightptr)
    {
      const Basetype *const midptr = leftptr
                                   + static_cast<size_t>(rightptr - leftptr)
                                                                / 2;
      if (pos < *midptr)
      {
        rightptr = midptr - 1;
      } else
      {
        found = midptr;
        if (pos > *midptr)
        {
          leftptr = midptr + 1;
        } else
        {
          return static_cast<size_t> (found - leftorig);
        }
      }
    }
    assert(found != nullptr && found >= leftorig);
    return size_t(1) + static_cast<size_t> (found - leftorig);
  }
 public:
  [[nodiscard]] size_t section_number_get(size_t elem) const noexcept
  {
    return elem >> logsectionsize;
  }
  [[nodiscard]] size_t section_min_elem(size_t idx) const noexcept
  {
    return idx << logsectionsize;
  }
  OrderedIntegerSequence(size_t _maxelement, size_t _nofelements)
    : logsectionsize(sizeof(Basetype) * CHAR_BIT)
    , nextfree(0)
    , maxelement(_maxelement)
    , currentsectionnum(0)
    , numofsections(section_number_get(_maxelement) + 1)
    , nofelements(_nofelements)
    , previouselem(undefined_element)
  {
    elements = new Basetype [nofelements];
    sectionstart = new size_t [numofsections + 1];
    sectionstart[0] = 0;
    for (size_t idx = 1; idx <= numofsections; idx++)
    {
      sectionstart[idx] = nofelements;
    }
  }
  ~OrderedIntegerSequence(void)
  {
    delete[] elements;
    delete[] sectionstart;
  }

  void show_info(void) const noexcept
  {
    std::cout << "# OrderedIntegerSequence(Basetype=";
    if (typeid(Basetype) == typeid(uint8_t))
    {
      std::cout << "uint8_t, ";
    } else
    {
      if (typeid(Basetype) == typeid(uint16_t))
      {
        std::cout << "uint16_t, ";
      } else
      {
        if (typeid(Basetype) == typeid(uint32_t))
        {
          std::cout << "uint32_t, ";
        }
      }
    }
    std::cout << "elements=" << nofelements
              << ", numofsections=" << numofsections << ", size="
              << (nofelements * sizeof(Basetype)
                  + (numofsections + 1) * sizeof(size_t))
              << " bytes)\n";
  }

  void append(size_t elem)
  {
    assert(nextfree < nofelements && elem <= maxelement &&
           (previouselem == undefined_element || previouselem < elem));

    while (elem >= section_min_elem(currentsectionnum + 1))
    {
      assert(currentsectionnum < numofsections);
      sectionstart[currentsectionnum + 1] = nextfree;
      currentsectionnum++;
    }
    assert(section_min_elem(currentsectionnum) <= elem &&
           elem < section_min_elem(currentsectionnum + 1) &&
           section_number_get(elem) == currentsectionnum);
    elements[nextfree++] = static_cast<Basetype>(elem);
    previouselem = elem;
  }

  [[nodiscard]] bool is_member(size_t elem) const noexcept
  {
    if (elem <= maxelement)
    {
      const size_t sectionnum = section_number_get(elem);

      if (sectionstart[sectionnum] < sectionstart[sectionnum + 1])
      {
        return binarysearch_is_member(elements + sectionstart[sectionnum],
                                      elements + sectionstart[sectionnum + 1]
                                               - 1,
                                      static_cast<Basetype>(elem));
      }
    }
    return false;
  }

  [[nodiscard]] size_t pos2seqnum(size_t pos) const noexcept
  {
    const size_t sectionnum = section_number_get(pos);
    assert(pos <= maxelement);

    if (sectionstart[sectionnum] < sectionstart[sectionnum + 1])
    {
      return sectionstart[sectionnum] +
             binarysearch_pos2seqnum(elements + sectionstart[sectionnum],
                                     elements + sectionstart[sectionnum+1] - 1,
                                     static_cast<Basetype>(pos));
    }
    return sectionstart[sectionnum];
  }

  void show_tabular(void) const noexcept
  {
    size_t sectionnum = 0;

    assert(nextfree > 0);
    for (size_t idx = 0; idx < nextfree; idx++)
    {
      while (idx >= sectionstart[sectionnum + 1])
      {
        sectionnum++;
      }
      std::cout << section_min_elem(sectionnum) + elements[idx]
                << (idx < nextfree - 1 ? "&" : "\\\\\n");
    }
    for (size_t idx = 0; idx < nextfree; idx++)
    {
      std::cout << elements[idx] << (idx < nextfree - 1 ? "&" : "\\\\\n");
    }
    sectionnum = 0;
    for (size_t idx = 0; idx < nextfree; idx++)
    {
      while (idx >= sectionstart[sectionnum + 1])
      {
        sectionnum++;
      }
      std::cout << sectionnum << (idx < nextfree - 1 ? "&" : "\\\\\n");
    }
    for (size_t idx = 0; idx <= numofsections; idx++)
    {
      std::cout << sectionstart[idx] << (idx < numofsections ? "&" : "\\\\\n");
    }
  }

  void show_text(void) const noexcept
  {
    assert(nextfree > 0);
    for (size_t sectionnum = 0; sectionnum < numofsections; sectionnum++)
    {
      std::cout << "# section\t" << sectionnum << "\t"
                << sectionstart[sectionnum] << '\n';
    }
    size_t sectionnum = 0;
    for (size_t idx = 0; idx < nextfree; idx++)
    {
      while (idx >= sectionstart[sectionnum + 1])
      {
        sectionnum++;
      }
      std::cout << "# elem\t" << idx << "\t"
                << sectionnum << "\t"
                << static_cast<unsigned int>(elements[idx]) << "\t"
                << (section_min_elem(sectionnum) + elements[idx])
                << std::endl;
    }
  }

  [[nodiscard]] size_t
  get_element_at(size_t sectionnum, size_t idx) const noexcept
  {
    assert(idx < nextfree);
    return section_min_elem(sectionnum) + elements[idx];
  }

  [[nodiscard]] size_t get_element_at(size_t idx) const noexcept
  {
    assert(idx < nextfree);
    const size_t sectionnum = binarysearch_sec_idx_largest_leq(idx);
#ifndef NDEBUG
    const size_t lsectionnum = linearsearch_sec_idx_largest_leq(idx);
    assert(sectionnum == lsectionnum);
#endif
    return section_min_elem(sectionnum) + elements[idx];
  }
};

size_t ordered_integer_sequence_size(size_t sizeofBasetype,
                                     size_t maxelement,
                                     size_t nofelements)
{
  const int logsectionsize = static_cast<int>(sizeofBasetype * CHAR_BIT);
  return nofelements * sizeofBasetype +
         ((maxelement >> logsectionsize) + 1) * sizeof(size_t);
}

size_t ordered_integer_sequence_sizeof_smallest(size_t maxelement,
                                                size_t nofelements)
{
  size_t sizeof_basetype = 0;
  size_t smallest_size = ~size_t(0);
  for (size_t this_sizeof = 1; this_sizeof <= 4; this_sizeof *= 2)
  {
    const size_t this_size
      = ordered_integer_sequence_size(this_sizeof,maxelement,nofelements);
    if (this_size < smallest_size)
    {
      sizeof_basetype = this_sizeof;
      smallest_size = this_size;
    }
  }
  assert(sizeof_basetype == 1 || sizeof_basetype == 2 || sizeof_basetype == 4);
  return sizeof_basetype;
}
#endif

#ifndef SORTED_MATCH_LIST_PP_HPP
#define SORTED_MATCH_LIST_PP_HPP

template<typename T>
struct SortedMatchListPositionPair
{
  T seqnum0,
    seqnum1,
    startpos1,
    startpos0;
  SortedMatchListPositionPair() {}
  SortedMatchListPositionPair(T _seqnum0,
                              T _seqnum1,
                              T _startpos1,
                              T _startpos0)
   : seqnum0(_seqnum0)
   , seqnum1(_seqnum1)
   , startpos1(_startpos1)
   , startpos0(_startpos0)
  {}
};
#endif

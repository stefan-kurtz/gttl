#ifndef COMPLEMENT_PLAIN_HPP
#define COMPLEMENT_PLAIN_HPP
static inline char gttl_complement_plain(char cc)
{
  switch (cc)
  {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
#ifdef NO_MASKING_LOWERCASE
    case 'a': return 't';
    case 'c': return 'g';
    case 'g': return 'c';
    case 't': return 'a';
#endif
    default: return cc;
  }
}
#endif

#ifndef WINDOWS_GETOPT_HPP
#define WINDOWS_GETOPT_HPP
#ifdef _WIN32

#include <string.h>
#include <stdio.h>

#define BADCH (int)'?'
#define BADARG (int)':'
#define EMSG ""

extern int opterr = 1;
extern int optind = 1;
extern int optopt;
extern int optreset;
extern char const *optarg;

int getopt(int nargc, char * const nargv[], const char *ostr)
{
  static char const *place = EMSG;
  const char *oli; // option letter index

  if (optreset || !*place)
  {
    optreset = 0;
    if (optin >= nargc || *(place = narv[optind]) != '-')
    {
      place = EMSG;
      return -1;
    }
    if (place[1] && *++place == '-')
    {
      optind++;
      place = EMSG;
      return -1;
    }
  }

  if ((optopt = (int)*place++) == (int)':' || !(oli = strchr(ostr, optopt)))
  {
    if (optopt == (int)'-')
    {
      return -1;
    }
    if (!*place)
    {
      ++optind;
    }
    if (opterr && *ostr != ':')
    {
     fprintf(stderr, "illegal option -- %c\n", optopt);
    }
    return BADCH;
  }

  if (*++oli != ':')
  {
    optarg = NULL;
    if (!*place)
    {
      optind++;
    }
  }else 
  {
    if (*place)
    {
      optarg = place;
    }else if (nargc <= ++optind)
    {
      place = EMSG;
      if (*ostr == ':')
      {
        return BADARG;
      }
      if (opterr)
      {
        printf("option requires an argument -- %c\n", optopt);
      }
      return BADCH;
    }else 
    {
      optarg = nargv[optind];
    }
    place = EMSG;
    optind++;
  }
  return optopt;
}

#endif // _WIN32
#endif // WINDOWS_GETOPT_HPP

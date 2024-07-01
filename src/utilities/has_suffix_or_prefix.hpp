#ifndef HAS_SUFFIX_OR_PREFIX_HPP
#define HAS_SUFFIX_OR_PREFIX_HPP
#include <string>
#include <vector>
static bool inline gttl_has_prefix(const std::string &s,
                                   const std::string &prefix)
{
  return s.size() >= prefix.size() and s.substr(0,prefix.size()) == prefix;
}

static bool inline gttl_has_suffix(const std::string &s,
                                   const std::string &suffix)
{
  return s.size() >= suffix.size() and
         s.substr(s.size() - suffix.size()) == suffix;
}

static inline bool gttl_likely_fasta_format(const std::string &inputfilename)
{
  const std::vector<std::string> fq_suffixes{".fq",".fastq",
                                             ".fq.gz",".fastq.gz"};
  for (auto &suffix : fq_suffixes)
  {
    if (gttl_has_suffix(inputfilename,suffix))
    {
      return false;
    }
  }
  return true;
}
#endif

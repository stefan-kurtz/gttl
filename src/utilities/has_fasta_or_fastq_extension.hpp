#ifndef HAS_FASTA_OR_FASTQ_EXTENSION_HPP
#define HAS_FASTA_OR_FASTQ_EXTENSION_HPP
#include <string>
#include <vector>
#include "utilities/has_suffix_or_prefix.hpp"

static inline bool gttl_likely_fasta_format(const std::string &inputfilename)
{
  const std::vector<std::string> fq_suffixes{".fq",".fastq",
                                             ".fq.gz",".fastq.gz"};

  for (const std::string& suf : fq_suffixes)
  {
    if (inputfilename.ends_with(suf)) return false;
  }
  return true;
}

static inline bool gttl_likely_gzipped_fastq_format(
  const std::string &inputfilename)
{
  return gttl_has_suffix_with_extension(inputfilename, ".fastq", ".gz") or
         gttl_has_suffix_with_extension(inputfilename, ".fq", ".gz");
}

static inline bool gttl_likely_gzipped_fasta_format(
  const std::string &inputfilename)
{
  // extensions according to https://en.wikipedia.org/wiki/FASTA_format
  static const std::vector<std::string> fasta_suffixes{".fasta",
                                                       ".fas",
                                                       ".fa",
                                                       ".fna",
                                                       ".ffn",
                                                       ".faa",
                                                       ".mpfa",
                                                       ".frn"};
  return gttl_has_any_suffix_with_extension(inputfilename,fasta_suffixes,".gz");
}
#endif

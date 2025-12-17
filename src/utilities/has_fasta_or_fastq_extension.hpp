#ifndef HAS_FASTA_OR_FASTQ_EXTENSION_HPP
#define HAS_FASTA_OR_FASTQ_EXTENSION_HPP
#include <algorithm>
#include <string>
#include <vector>
#include "utilities/has_suffix_or_prefix.hpp"

static inline bool gttl_likely_fasta_format(const std::string &inputfilename)
{
  const std::vector<std::string> fq_suffixes{".fq",".fastq",
                                             ".fq.gz",".fastq.gz"};
  return std::ranges::none_of(fq_suffixes,
                              [&inputfilename](const std::string& suf)
                              {
                                return inputfilename.ends_with(suf);
                              });
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

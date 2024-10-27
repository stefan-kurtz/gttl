#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <zlib.h>
#include "sequences/gttl_fastq_iterator.hpp"
#include "sequences/gttl_seq_iterator.hpp"
#include "sequences/split_files.hpp"
#include "utilities/cxxopts.hpp"
#include "utilities/gttl_line_iterator.hpp"
#include "utilities/has_suffix_or_prefix.hpp"

int main(int argc, char *argv[])
{
  cxxopts::Options options("split_gzip_mn.x",
                           "Split gzip-compressed sequence files");
  options.set_width(80);
  options.set_tab_expansion();

  options.add_options()("f,file", "Input file to split",
                        cxxopts::value<std::string>())(
      "h,help", "Print usage information")(
      "n,num_parts",
      "Number of parts to split the file into (mutually exclusive with -l, -s)",
      cxxopts::value<size_t>())(
      "l,len_parts",
      "Length of parts to split the file into (mutually exclusive with -n, -s)",
      cxxopts::value<size_t>())(
      "s,num_sequences",
      "Number of sequences per file (mutually exclusive with -l, -n)",
      cxxopts::value<size_t>())(
      "o,output",
      "The output filenames (numbering and file-extensions will be added)",
      cxxopts::value<std::string>())(
      "c,compression_level", "The GZip compression level of the output files",
      cxxopts::value<size_t>()->default_value("6"));

  options.parse_positional({"file"});
  options.positional_help("<input file>");

  auto result = options.parse(argc, argv);

  if (result.count("help"))
  {
    std::cout << options.help() << std::endl;
    exit(EXIT_SUCCESS);
  }

  if (!result.count("file"))
  {
    std::cout << options.help() << std::endl;
    exit(EXIT_FAILURE);
  }

  if (result.count("num_parts") + result.count("len_parts") +
          result.count("num_sequences") !=
      1)
  {
    std::cerr << argv[0] << ": Exactly one of -l and -n must be specified"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  size_t num_parts, len_parts, num_sequences;
  if (result.count("num_parts"))
  {
    num_parts = result["num_parts"].as<size_t>();
    len_parts = 0;
    num_sequences = 0;
  } else if (result.count("len_parts"))
  {
    num_parts = 0;
    len_parts = result["len_parts"].as<size_t>();
    num_sequences = 0;
  } else
  {
    num_parts = 0;
    len_parts = 0;
    num_sequences = result["num_sequences"].as<size_t>();
  }

  std::string ifilename = result["file"].as<std::string>();

  std::string output_basename;
  if (result.count("output"))
  {
    output_basename = result["output"].as<std::string>();
  } else
  {
    output_basename = ifilename.substr(ifilename.find_last_of("/") + 1,
                                       ifilename.find_last_of("."));
    if (gttl_has_suffix(ifilename,
                        ".gz"))  // We need to remove two file extensions here
    {
      output_basename =
          output_basename.substr(0, output_basename.find_last_of("."));
    }
    output_basename += "_part";
  }

  constexpr const int BUF_SIZE = 1 << 14;

  if (gttl_has_suffix(ifilename, ".fasta") ||
      gttl_has_suffix(ifilename, ".fasta.gz"))
  {
    GttlSeqIterator<BUF_SIZE> fasta_it(ifilename.c_str());
    if (num_parts != 0)
    {
      split_into_num_files(fasta_it, output_basename, num_parts,
                           result["compression_level"].as<size_t>());
      return EXIT_SUCCESS;
    }
    if (len_parts != 0)
    {
      split_into_parts_length(fasta_it, output_basename, len_parts,
                              result["compresssion_level"].as<size_t>());
      return EXIT_SUCCESS;
    }
    if (num_sequences != 0)
    {
      split_into_num_sequences(fasta_it, output_basename, num_sequences,
                               result["compression_level"].as<size_t>());
      return EXIT_SUCCESS;
    }
  } else if (gttl_has_suffix(ifilename, ".fastq") ||
             gttl_has_suffix(ifilename, ".fastq.gz"))
  {
    GttlLineIterator<BUF_SIZE> line_it(ifilename.c_str());
    GttlFastQIterator<GttlLineIterator<BUF_SIZE>> fastq_it(line_it);
    if (num_parts != 0)
    {
      split_into_num_files(fastq_it, output_basename, num_parts,
                           result["compression_level"].as<size_t>());
      return EXIT_SUCCESS;
    }
    if (len_parts != 0)
    {
      split_into_parts_length(fastq_it, output_basename, len_parts,
                              result["compression_level"].as<size_t>());
      return EXIT_SUCCESS;
    }
    if (num_sequences != 0)
    {
      split_into_num_sequences(fastq_it, output_basename, num_sequences,
                               result["compression_level"].as<size_t>());
      return EXIT_SUCCESS;
    }
  }

  std::cerr << argv[0] << ": file " << ifilename
            << " does not have .fasta[.gz]|.fastq[.gz] extension!" << std::endl;
  exit(EXIT_FAILURE);
}

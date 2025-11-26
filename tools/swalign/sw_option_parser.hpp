#ifndef SW_OPTION_PARSER_HPP
#define SW_OPTION_PARSER_HPP
#ifndef _WIN32
  // We disable the include-cleaner check here, since
  // optind is defind in <bits/getopt_core.h>, which should *NOT* be
  // included directly.
  #include <unistd.h> // NOLINT(misc-include-cleaner)
#else
  #define NOMINMAX
  #include <io.h>
  #include "utilities/windows_getopt.hpp"
#endif
#include <cstdint>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <exception>
#include <stdexcept>
#include <format>
#include "threading/threads_output_files.hpp"
#include "alignment/score_matrix_name.hpp"
#include "alignment_display.hpp"

struct SWOptions
{
  static constexpr const size_t min_alignment_width = size_t(8);
  bool sw_all_against_all;
  int8_t gap_open_penalty,
         gap_extension_penalty;
  char *threads_out_prefix, *restrict_to_pairs_file, *dbfile, *queryfile;
  bool header_display,
       opt_memory,
       stop_after_first,
       no_reverse_strand,
       all_against_all;
  int vectorized_alignment; /* 0: none
                               1: raw_score and end-positions
                               2: raw_score, end-positions, start-positions */
  AlignmentDisplay alignment_display;
  size_t best, num_threads;
  uint32_t min_bit_score;
  const char *score_matrix_id;
  ScoreMatrixName score_matrix_name;
  SWOptions(bool _sw_all_against_all)
    : sw_all_against_all(_sw_all_against_all)
    , gap_open_penalty(int8_t(11))
    , gap_extension_penalty(int8_t(1))
    , threads_out_prefix(nullptr)
    , restrict_to_pairs_file(nullptr)
    , dbfile(nullptr)
    , queryfile(nullptr)
    , header_display(false)
    , opt_memory(false)
    , stop_after_first(false)
    , no_reverse_strand(false)
    , all_against_all(true)
    , vectorized_alignment(2)
    , alignment_display({})
    , best(0)
    , num_threads(size_t(1))
    , min_bit_score(0)
    , score_matrix_id(nullptr)
    , score_matrix_name({})
  {}

  void parse(int argc, char * const argv[])
  {
    assert(argc > 1);
    while (true)
    {
      const char *const flags = sw_all_against_all ? "dqsgvcthnbamfor" :
                                                     "dqsgvcthnbamp";
      const int opt = getopt(argc,argv,flags);
      const int lower_bound_vectorized = sw_all_against_all ? 1 : 0;
      if (opt == -1)
      {
        break;
      }
      const char c_opt = static_cast<char>(opt);
      switch (c_opt)
      {
        case 'd':
          if (optind > argc - 1)
          {
            throw std::invalid_argument("missing argument to option -d");
          }
          dbfile = argv[optind++];
          break;
        case 'q':
          if (optind > argc - 1)
          {
            throw std::invalid_argument("missing argument to option -q");
          }
          queryfile = argv[optind++];
          break;
        case 'g':
          {
            int readint;
            if (optind + 1 > argc - 1)
            {
              throw std::invalid_argument("missing arguments to option -g");
            }
            if (std::sscanf(argv[optind],"%d",&readint) != 1 || readint < 0 ||
                            readint > INT8_MAX)
            {
              throw std::invalid_argument(
                "illegal first argument to option -g");
            }
            gap_open_penalty = static_cast<int8_t>(readint);
            if (std::sscanf(argv[optind+1],"%d",&readint) != 1 || readint < 0 ||
                            readint > INT8_MAX)
            {
              throw std::invalid_argument(
                "illegal second argument to option -g");
            }
            gap_extension_penalty = static_cast<int8_t>(readint);
            optind += 2;
            break;
          }
        case 't':
        case 'b':
        case 'c':
          {
            int readint;
            if (optind > argc - 1 ||
                std::sscanf(argv[optind],"%d",&readint) != 1)
            {
              throw std::invalid_argument(
                      std::format("missing or illegal argument to option -{}",
                                  c_opt));
            }
            if (c_opt == 't')
            {
              if (readint < 1)
              {
                throw std::invalid_argument(
                  "argument to option -t (i.e. the number of "
                  "threads) must be positive");
              }
              num_threads = static_cast<size_t>(readint);
            } else
            {
              if (c_opt == 'b')
              {
                if (readint < 1)
                {
                  throw std::invalid_argument(
                    "argument to option -b must be positive");
                }
                best = static_cast<size_t>(readint);
              } else
              {
                if (readint < 1)
                {
                  throw std::invalid_argument(
                    "argument to option -c must be positive");
                }
                assert(c_opt == 'c');
                min_bit_score = static_cast<uint32_t>(readint);
              }
            }
            optind++;
            break;
          }
        case 'a':
          {
            if (optind <= argc - 1)
            {
              if (strlen(argv[optind]) > 0 and argv[optind][0] != '-')
              {
                alignment_display.set_from_string(min_alignment_width,
                                                  argv[optind]);
              } else
              {
                if (strlen(argv[optind]) > 0 and argv[optind][0] == '-')
                {
                  throw std::invalid_argument("missing argument to option -a");
                }
              }
              optind++;
            } else
            {
              throw std::invalid_argument("missing argument to option -a");
            }
            break;
          }
        case 's':
        case 'o':
        case 'r':
          {
            if (optind > argc - 1)
            {
              throw std::invalid_argument(
                      std::format("missing argument to option -{}", c_opt));
            }
            if (c_opt == 's')
            {
              score_matrix_id = argv[optind++];
              try
              {
                score_matrix_name.set(std::string(score_matrix_id));
              }
              catch(const std::exception &err)
              {
                throw;
              }
            } else
            {
              if (c_opt == 'o')
              {
                threads_out_prefix = argv[optind++];
              } else
              {
                assert(c_opt == 'r');
                restrict_to_pairs_file = argv[optind++];
              }
            }
            break;
          }
        case 'v':
          {
            int readint;
            if (optind > argc - 1 ||
                std::sscanf(argv[optind],"%d",&readint) != 1 ||
                readint < lower_bound_vectorized ||
                readint > 2)
            {
              std::string msg{"missing or illegal argument to option -v; "
                              "argument must be either "};
              msg += sw_all_against_all ? "1 or 2" :  "0, 1, or 2";
              throw std::invalid_argument(msg);
            }
            vectorized_alignment = readint;
            optind++;
            break;
          }
        case 'h':
          header_display = true;
          break;
        case 'n':
          no_reverse_strand = true;
          break;
        case 'p':
          all_against_all = false;
          break;
        case 'm':
          opt_memory = true;
          break;
        case 'f':
          stop_after_first = true;
          break;
        default:
          {
          throw std::invalid_argument(std::format("illegal option -{}", c_opt));
          }
      }
    }
    if (optind < argc)
    {
      throw std::invalid_argument(std::format("superfluous arguments {}",
                                              argv[optind]));
    }
    assert (optind == argc);
    if (dbfile == nullptr)
    {
      throw std::invalid_argument("option -d is mandatory");
    }
    if (vectorized_alignment == 1 && alignment_display.need_alignment())
    {
      vectorized_alignment = 2;
    }
    if (!all_against_all)
    {
      if (best > 0)
      {
        throw std::invalid_argument("option -p and -b are incompatible");
      }
      if (queryfile != nullptr)
      {
        throw std::invalid_argument("option -p and -q are incompatible");
      }
      if (num_threads > 1)
      {
        throw std::invalid_argument("option -p and -t are incompatible");
      }
    }
    if (queryfile == nullptr)
    {
      queryfile = dbfile;
    }
  }

  void show(FILE *fpout,const char *progname) const noexcept
  {
    fprintf(fpout,"# Options: %s",progname);
    if (threads_out_prefix != nullptr)
    {
      fprintf(fpout," -o %s",threads_out_prefix);
    }
    if (restrict_to_pairs_file != nullptr)
    {
      fprintf(fpout," -r %s",restrict_to_pairs_file);
    }
    fprintf(fpout," -v %d",vectorized_alignment);
    fprintf(fpout," -g %d %d",static_cast<int>(gap_open_penalty),
                              static_cast<int>(gap_extension_penalty));
    if (alignment_display.need_alignment())
    {
      fprintf(fpout," -a %s",alignment_display.to_string().c_str());
    }
    if (best > 0)
    {
      fprintf(fpout," -b %zu",best);
    }
    if (min_bit_score > 0)
    {
      fprintf(fpout," -c %u",min_bit_score);
    }
    if (header_display)
    {
      fprintf(fpout," -h");
    }
    if (num_threads > 1)
    {
      fprintf(fpout," -t %zu",num_threads);
    }
    if (dbfile != nullptr)
    {
      fprintf(fpout," -d %s",dbfile);
    }
    if (queryfile != nullptr && strcmp(dbfile,queryfile) != 0)
    {
      fprintf(fpout," -q %s",queryfile);
    }
    fprintf(fpout,"\n");
  }
  void usage(const char *progname) const noexcept
  {
    fprintf(stderr,"Usage: %s [options]\n",progname);
    fprintf(stderr,"options:\n");
    if (not sw_all_against_all)
    {
      fprintf(stderr,
    "    -p compare pairs of consecutive sequences rather than all\n"
    "       against all, which is the default comparison strategy;\n"
    "       this option is incompatible with option -b, -t, -q\n");
    }
    fprintf(stderr,
    "    -d <dbfile> (mandatory)\n"
    "    -q <queryfile> (optional)\n"
    "    -s <scorematrix> (default: blosum62 for proteins, unit_score_nuc for\n"
    "                      DNA sequences, other matrices for proteins have\n"
    "                      suffix _aa; matrices for DNA sequences have suffix\n"
    "                      _nuc; possible values: %s)\n",
        score_matrix_name.string_values_joined(", ").c_str());
    fprintf(stderr,
    "    -g <gap_open_penalty> <gap_extension_penalty> (default: 11 1)\n"
    "    -v mode of computation\n");
    if (not sw_all_against_all)
    {
      fprintf(stderr,
      "          0: full dynamic programming, not vectorized and thus slow\n"
      "             (use for test purpose only)\n");
    }
    fprintf(stderr,
    "          1: vectorized, reporting score and end positions\n"
    "          2: vectorized, reporting score and start + end positions\n"
    "             (default)\n");
    fprintf(stderr,
    "    -n no_reverse_strand: when processing DNA sequences do not perform\n"
    "       comparison, in which one sequence is from the reverse strand\n"
    "       (default: false)\n");
    fprintf(stderr,
    "    -h switch to display sequence header (default: false)\n"
    "    -c <minimum bit score of alignments displayed,\n"
    "       if bit score is not available filter for raw score> (default: 0)\n"
    "    -t <number of threads>\n");
    if (sw_all_against_all)
    {
      fprintf(stderr,"    -o <thread_out_prefix> (optional); "
                     "%s (default \"\")\n",
              ThreadsOutputFiles::help_line);
    }
    fprintf(stderr,"    -r <inputfile> "
                   "specify file with pairs of sequence headers\n"
    "                   (one per line) to restrict computation to\n");
    fprintf(stderr,
    "    -a <number> compute and possibly display alignment\n"
    "          1: verify score of aligned regions (using banded DP-method)\n"
    "          2: display coverage of aligned subject sequence relative to\n"
    "             complete subject sequence,\n"
    "          3: display coverage of aligned query sequence relative to\n"
    "             complete query sequence,\n"
    "          4: compute alignment using banded method and display\n"
    "             identity,\n"
    "          5: compute alignment using banded method and display\n"
    "             corresponding cigar string,\n"
    "          6: output subject sequence involved in alignment,\n"
    "          7: output query sequence involved in alignment,\n"
    "          >=%zu: compute alignment and display it in C columns\n"
    "          where C is the number given as argument.\n"
    "          One can concatenate ordered values\n"
    "          >= 2 with a +-sign to specify their combination.\n"
    "          For example, -a 4+60 means to show the\n"
    "          identity and the alignment in 60 columns; -a 4+5+6 means to\n"
    "          show the identity, the cigar string and the subject sequences;\n"
    "          default: 0\n",
            min_alignment_width);
    fprintf(stderr,
    "    -m trade less space for more time when displaying alignments\n"
    "       (default: false)\n"
    "    -b <number of best local alignments to show> (if not used, then for\n"
    "       all pairs of sequences an optimal local alignment is shown)\n"
    );
  }
  void show_fields(FILE *fpout,
                   bool dna_alphabet,
                   bool with_bit_score)
  {
    const char *const tag = header_display ? "id" : "num";

    if (vectorized_alignment == 0 || vectorized_alignment == 2)
    {
      fprintf(fpout,"# Fields: s. seq%s, q. seq%s, s. start, s. len, "
                    "q. start, q. len, score",tag,tag);
    } else
    {
      assert(vectorized_alignment == 1);
      fprintf(fpout,"# Fields: s. seq%s, q. seq%s, s. end, q. end, score",
              tag,tag);
    }
    if (with_bit_score > 0)
    {
      fprintf(fpout,", bit score");
    }
    if (dna_alphabet)
    {
      fprintf(fpout,", strand");
    }
    if (alignment_display.s_coverage())
    {
      fprintf(fpout,", s. cov");
    }
    if (alignment_display.q_coverage())
    {
      fprintf(fpout,", q. cov");
    }
    if (alignment_display.identity())
    {
      fprintf(fpout,", identity");
    }
    if (alignment_display.cigar())
    {
      fprintf(fpout,", cigar");
    }
    if (alignment_display.s_substring())
    {
      fprintf(fpout,", s. substr");
    }
    if (alignment_display.q_substring())
    {
      fprintf(fpout,", q. substr");
    }
    fputc('\n',fpout);
  }
};
#endif

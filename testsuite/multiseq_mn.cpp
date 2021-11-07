#include <getopt.h>
#include <iostream>
#include <algorithm>
#include <array>
#include "utilities/runtime_class.hpp"
#include "sequences/gttl_multiseq.hpp"
#include "sequences/literate_multiseq.hpp"

static void usage_output(bool error_case,const char *progname)
{
  StrFormat msg("Usage: %s [options] <filename1> [filename2 ..]\n"
                "  -h --help            Show this help\n"
                "  -p --protein         handle protein sequences\n"
                "  -w --width <width>   output headers and sequences; \n"
                "                       width specifies the line width of the\n"
                "                       sequence output; 0 means to output\n"
                "                       sequence in a single line",progname);
  if (error_case)
  {
    std::cerr << msg.str() << std::endl;
  } else
  {
    std::cout << msg.str() << std::endl;
  }
}

int main(int argc, char *argv[])
{
  /* Different variables used for the optionparser as well as multiseq variable
     multiseq */
  /* Give help when no arguments were given */
  if (argc == 1)
  {
    usage_output(true,argv[0]);
    return EXIT_FAILURE;
  }

  /* getopt commandline parser */
  const struct option longopts[] = {
    {"help", no_argument, 0, 'h'},
    {"width", required_argument, 0, 'w'},
    {"protein", no_argument, 0, 'p'},
    {"rankdist", no_argument, 0, 'r'},
    {0, 0, 0, 0}
  };

  int width = -1;
  bool protein = false, rankdist = false;
  int opt;
  while ((opt = getopt_long(argc, argv, ":hprw:", longopts, NULL)) != -1)
  {
    switch (opt)
    {
      case 'h':
        usage_output(false,argv[0]);
        return EXIT_SUCCESS;
        break;
      case 'p':
        protein = true;
        break;
      case 'r':
        rankdist = true;
        break;
      case 'w':
        if (sscanf(optarg,"%d",&width) != 1 || width < 0)
        {
          std::cerr << argv[0] << ": argument to option -"
                    << static_cast<char>(opt)
                    << " must be non-negative integer"
                    << std::endl;
          return EXIT_FAILURE;
        }
        break;
      case 0:
        break;
      case ':':
        std::cerr << argv[0] << ": Option -" << optopt
                             << " requires an argument" << std::endl;
        usage_output(true,argv[0]);
        return EXIT_FAILURE;
    }
  }
  if (optind > argc-1)
  {
    std::cerr << argv[0] << ": missing filename argument" << std::endl;
    usage_output(true,argv[0]);
    return EXIT_FAILURE;
  }
  GttlMultiseq *multiseq = nullptr;
  RunTimeClass rt = RunTimeClass();
  std::vector<std::string> inputfiles{};
  try
  {
    for (int idx = optind; idx < argc; idx++)
    {
      inputfiles.push_back(std::string(argv[idx]));
    }
    const bool store_sequences = (width >= 0 || rankdist) ? true : false;
    multiseq = new GttlMultiseq(inputfiles,store_sequences,UINT8_MAX);
  }
  catch (std::string &msg)
  {
    for (auto &&inputfile : inputfiles)
    {
      std::cerr << argv[0] << ": file \"" << inputfile << "\""
                << msg << std::endl;
    }
    delete multiseq;
    return EXIT_FAILURE;
  }
  rt.show("create GttlMultiseq");
  if (width >= 0)
  {
    multiseq->show(static_cast<size_t>(width), false);
  }
  for (auto &&inputfile : inputfiles)
  {
    std::cout << "# filename\t" << inputfile << std::endl;
  }
  std::cout << "# number of sequences\t"
            << multiseq->sequences_number_get()
            << std::endl;
  std::cout << "# total length\t"
            << multiseq->sequences_total_length_get()
            << std::endl;
  if (rankdist)
  {
    if (protein)
    {
      static constexpr const char amino_acids[]
        = "A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y";
      LiterateMultiseq<amino_acids,20> literate_multiseq(*multiseq);
      literate_multiseq.show_rank_dist();
    } else
    {
      static constexpr const char nucleotides_upper_lower[] = "Aa|Cc|Gg|TtUu";
      LiterateMultiseq<nucleotides_upper_lower,4> literate_multiseq(*multiseq);
      literate_multiseq.show_rank_dist();
    }
  }
  delete multiseq;
}

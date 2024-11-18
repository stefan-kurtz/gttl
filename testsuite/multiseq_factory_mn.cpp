#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <cinttypes>
#include "utilities/runtime_class.hpp"
#include "sequences/multiseq_factory.hpp"

int main(int argc, char *argv[])
{
  int64_t number_of_sequences_in_split;
  GttlMultiseqFactory *multiseq_factory = nullptr;

  if ((argc == 3 or argc == 4) and
      std::sscanf(argv[1],"%" PRId64,&number_of_sequences_in_split) == 1 and
      number_of_sequences_in_split > 0)
  {
    try
    {
      const bool short_header = true;
      const uint8_t padding_char = UINT8_MAX;
      if (argc == 3)
      {
        multiseq_factory
          = new GttlMultiseqFactory(std::string(argv[2]),
                                    static_cast<size_t>
                                               (number_of_sequences_in_split),
                                    padding_char,
                                    short_header);
      } else
      {
        multiseq_factory
          = new GttlMultiseqFactory(std::string(argv[2]),
                                    std::string(argv[3]),
                                    static_cast<size_t>
                                               (number_of_sequences_in_split),
                                    padding_char,
                                    short_header);
      }
    }
    catch (std::string &msg)
    {
      std::cerr << argv[0] << ": file \"" << argv[2] << "\"";
      if (argc == 4)
      {
        std::cerr << ", \"" << argv[2] << "\"";
      }
      std::cerr << msg << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << "# number of parts\t" << multiseq_factory->size() << std::endl;
    delete multiseq_factory;
  } else
  {
    std::cerr << "Usage: " << argv[0]
              << " num_seq_in_split filename0 [filename1]" << std::endl;
    return EXIT_FAILURE;
  }
}

#include <iostream>
#include <vector>
#include <cinttypes>
#include "sequences/multiseq_factory.hpp"

template<bool length_based>
static void test_multiseq_factory(size_t number_of_sequences_in_split,
                                  const std::vector<std::string> &inputfiles)
{
  const uint8_t padding_char = UINT8_MAX;
  const bool short_header = true;
  GttlMultiseqFactory<length_based> *multiseq_factory = nullptr;
  if (inputfiles.size() == 1)
  {
    multiseq_factory
      = new GttlMultiseqFactory<length_based>
                               (inputfiles[0],
                                static_cast<size_t>
                                           (number_of_sequences_in_split),
                                padding_char,
                                short_header);
  } else
  {
    assert (inputfiles.size() == 2);
    multiseq_factory
      = new GttlMultiseqFactory<length_based>
                               (inputfiles[0],
                                inputfiles[1],
                                static_cast<size_t>
                                           (number_of_sequences_in_split),
                                padding_char,
                                short_header);
  }
  std::cout << "# number of parts\t" << multiseq_factory->size() << std::endl;
  delete multiseq_factory;
}

int main(int argc, char *argv[])
{
  int64_t number_of_sequences_in_split;
  static constexpr const bool length_based = false;

  if ((argc == 3 or argc == 4) and
      std::sscanf(argv[1],"%" PRId64,&number_of_sequences_in_split) == 1 and
      number_of_sequences_in_split > 0)
  {
    try
    {
      std::vector<std::string> inputfiles{std::string(argv[2])};
      if (argc == 4)
      {
        inputfiles.push_back(std::string(argv[3]));
      }
      test_multiseq_factory<length_based>(number_of_sequences_in_split,
                                          inputfiles);
    }
    catch (std::string &msg)
    {
      std::cerr << argv[0] << ": file \"" << argv[2] << "\"";
      if (argc == 4)
      {
        std::cerr << ", \"" << argv[4] << "\"";
      }
      std::cerr << msg << std::endl;
      return EXIT_FAILURE;
    }
  } else
  {
    std::cerr << "Usage: " << argv[0]
              << " num_seq_in_split filename0 [filename1]" << std::endl;
    return EXIT_FAILURE;
  }
}

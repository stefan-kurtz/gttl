#include <cstdlib>
#include <iostream>
#include <cstdint>
#include "utilities/read_vector.hpp"

int main(int argc,char *argv[])
{
  if (argc != 3)
  {
    std::cout << "Usage: " << argv[0] << " file file.gz" << std::endl;
    std::cout << "such that file.gz is the compressed version of file"
              << std::endl;
    return EXIT_FAILURE;
  }
  try
  {
    const auto file_content = gttl_read_vector<uint8_t>(argv[1]),
               file_content_from_zipped = gttl_read_vector<uint8_t>(argv[2]);
    if (file_content.size() != file_content_from_zipped.size())
    {
      std::cerr << argv[0] << ": "
                << "file " << argv[1]
                << " (size " << file_content.size() << ") and "
                << "file " << argv[2]
                << " (size " << file_content_from_zipped.size()
                << ") differ in size" << std::endl;
      return EXIT_FAILURE;
    }
    if (file_content != file_content_from_zipped)
    {
      std::cerr << argv[0] << ": " << "file content of " << argv[1]
                << " and " << argv[2] << " differ" << std::endl;
      return EXIT_FAILURE;
    }
  } catch(const std::string &msg)
  {
    std::cerr << argv[0] << ": " << msg << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

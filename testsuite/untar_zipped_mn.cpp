#include <cstdlib>
#include <exception>
#include <iostream>
#include <cstdio>
#include "utilities/split_string.hpp"
#include "utilities/is_in_PATH.hpp"
#include "utilities/untar_zipped.hpp"
#include "untar_zipped_op.hpp"

static void entry_show(const DecompressedFile &entry)
{
  std::cout << entry.filename_get() << "\t"
            << entry.size() << "\t"
            << (entry.is_directory() ? "d" : "f")
            << std::endl;
  if (entry.size() < 1000)
  {
    std::cout << "'''file_contents" << std::endl;
    std::cout << std::string_view((const char *) entry.data(),entry.size());
    std::cout << "'''" << std::endl;
  }
}

int main(int argc, char* argv[])
{
  UnzippedTarOptions options;
  try
  {
    options.parse(argc,argv);
  }
  catch (const std::invalid_argument &e)
  {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  if (options.help_option_is_set())
  {
    return EXIT_SUCCESS;
  }
  std::vector<DecompressedFile> decompressed_files;
  try
  {
    const bool with_rapidgzip = (not options.no_rapidgzip_option_is_set()) and
                                gttl_is_in_PATH("gtar") and
                                gttl_is_in_PATH("rapidgzip");
    for (auto &&inputfile : options.inputfiles_get())
    {
      TarReader tar_reader(inputfile,with_rapidgzip,false);
      for (auto entry : tar_reader)
      {
        if (options.store_option_is_set())
        {
          decompressed_files.push_back(entry);
        } else
        {
          entry_show(entry);
        }
      }
    }
  }
  catch (const std::exception& err)
  {
    std::cerr << argv[0] << ": " << err.what() << std::endl;
    return EXIT_FAILURE;
  }
  if (options.store_option_is_set())
  {
    for (auto &entry: decompressed_files)
    {
      entry_show(entry);
    }
  }
  return EXIT_SUCCESS;
}

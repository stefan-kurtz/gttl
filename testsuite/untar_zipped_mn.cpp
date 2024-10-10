#include <cstdlib>
#include <iostream>
#include <cstdio>
#include "utilities/untar_zipped.hpp"

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
  if (argc != 3 or (std::string(argv[1]) != "store" and
                    std::string(argv[1]) != "stream"))
  {
    std::cerr << "Usage: " << argv[0]
              << " <store|stream> <inputile (.tar.gz/.tar.bz2)>"
              << std::endl;
    return EXIT_FAILURE;
  }
  const char *inputfile = argv[2];

  const bool stream = std::string(argv[1]) == "stream";
  std::vector<DecompressedFile> decompressed_files;
  try
  {
    TarReader tar_reader(inputfile);

    for (auto entry : tar_reader)
    {
      if (stream)
      {
        entry_show(entry);
        entry.delete_data();
      } else
      {
        decompressed_files.push_back(entry);
      }
    }
  }
  catch (const std::string& msg)
  {
    std::cerr << argv[0] << ": " << msg << std::endl;
    return EXIT_FAILURE;
  }
  if (not stream)
  {
    for (auto &entry: decompressed_files)
    {
      entry_show(entry);
      entry.delete_data();
    }
  }
  return EXIT_SUCCESS;
}

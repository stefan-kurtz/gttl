#include <algorithm>
#include <cstring>
#include <iostream>
#include "sequences/gttl_fastq_generator.hpp"
#include "utilities/gttl_file_open.hpp"

#define FILENAME "../testdata/SRR19536726_1_1000.fastq.gz"

inline size_t count_occ(const char* string, char symbol)
{
  return std::count(string, string+std::strlen(string), symbol);
}

int main()
{
  size_t pseudo_hash = 0;
  GttlFastQGenerator fg(gttl_fp_type_open(FILENAME, "rb"));
  for(auto entry : fg)
  {
    pseudo_hash += count_occ(entry->header, '@');
    pseudo_hash += count_occ(entry->sequence, 'a');
    pseudo_hash += count_occ(entry->quality, 'F');
  }
  std::cout << pseudo_hash << "\n";
}

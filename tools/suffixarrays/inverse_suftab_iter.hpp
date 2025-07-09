#ifndef INVERSE_SUFTAB_ITER_HPP
#define INVERSE_SUFTAB_ITER_HPP

#include <cstddef>
#include <cstdio>
#include <ios>
#include <string>
#include <cassert>
#include <filesystem>
#include "utilities/gttl_binary_read.hpp"
#include "utilities/memory_tracker.hpp"

template <typename SuftabBaseType>
class InverseSuftabReader
{
  const std::string inverse_suftab_outputfile;
  const BinaryFileReader<SuftabBaseType> inverse_suftab_reader;
  [[nodiscard]] bool isa_file_is_up_to_date(const std::string &indexname) const
  {
    const std::string suftab_inputfile(indexname + ".suf");
    const std::string inverse_suftab_inputfile(indexname + ".isa");
    return std::filesystem::exists(suftab_inputfile) and
           std::filesystem::exists(inverse_suftab_inputfile) and
           std::filesystem::last_write_time(suftab_inputfile) <
           std::filesystem::last_write_time(inverse_suftab_inputfile);
  }
  const std::string create_inverse_suftab_file(
              GttlMemoryTracker *memory_tracker,
              const std::string &indexname,
              size_t totallength) const
  {
    const std::string inverse_suftab_filename{indexname + ".isa"};
    if (not isa_file_is_up_to_date(indexname))
    {
      const std::string suftab_inputfile(indexname + ".suf");
      const BinaryFileReader<SuftabBaseType> suftab_reader(suftab_inputfile);
      auto suftab_iter = suftab_reader.begin();
      SuftabBaseType *const inverse_suftab =
                                   new SuftabBaseType[totallength + 1];
      memory_tracker->track(inverse_suftab, __FILE__, __LINE__,
                            (totallength + 1) * sizeof(SuftabBaseType));
      for (size_t idx = 0; idx <= totallength; idx++)
      {
        inverse_suftab[*suftab_iter] = idx;
        assert(suftab_iter != suftab_reader.end());
        ++suftab_iter;
      }
      assert(suftab_iter == suftab_reader.end());
      FILE *const out_fp = std::fopen(inverse_suftab_filename.c_str(), "wb");
      if (out_fp == nullptr)
      {
        throw std::ios_base::failure(std::string("cannot create file ") +
                                     inverse_suftab_filename);
      }
      std::fwrite(inverse_suftab, sizeof *inverse_suftab, totallength + 1,
                  out_fp);
      std::fclose(out_fp);
      memory_tracker->untrack(inverse_suftab, __FILE__, __LINE__);
      delete[] inverse_suftab;
    }
    return inverse_suftab_filename;
  }
  public:
  InverseSuftabReader(GttlMemoryTracker *memory_tracker,
                      const std::string &indexname, size_t totallength)
    : inverse_suftab_outputfile(create_inverse_suftab_file(memory_tracker,
                                                           indexname,
                                                           totallength))
    , inverse_suftab_reader(inverse_suftab_outputfile)
  { }
  [[nodiscard]] auto begin(void) const { return inverse_suftab_reader.begin(); }
  [[nodiscard]] auto end(void) const { return inverse_suftab_reader.end(); }
  ~InverseSuftabReader(void)
  { }
};
#endif

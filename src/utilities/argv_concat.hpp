#ifndef ARGV_CONCAT_HPP
#define ARGV_CONCAT_HPP
#include <string.h> // NOLINT(modernize-deprecated-headers)
#include <cstdlib>
#include <string>
#include <vector>
#include <cstring>
#include <cassert>
#include <iostream>

class ArgvConcat
{
  private:
  std::vector<std::string> arg_vector;
  char **my_argv;
  public:
  ArgvConcat(int argc,const char* const* argv,const char *separator,
             bool replace_single_hyphen_options,
             const std::vector<std::string> &concat_options) :
    arg_vector({}),
    my_argv(nullptr)
  {
    bool concat = false;
    for (int idx = 0; idx < argc; idx++)
    {
      if (concat)
      {
        assert(not arg_vector.empty());
        const std::string key(argv[idx]);
        if (strlen(argv[idx]) > 0 && argv[idx][0] != '-')
        {
          if (arg_vector.back().at(0) == '-')
          {
            arg_vector.emplace_back(argv[idx]);
          } else
          {
            arg_vector.back() += (separator + std::string(argv[idx]));
          }
        } else
        {
          concat = false;
          arg_vector.emplace_back(argv[idx]);
        }
      } else
      {
        arg_vector.emplace_back(argv[idx]);
        for (auto &&opt : concat_options)
        {
          if (std::string(argv[idx]) == opt)
          {
            concat = true;
            break;
          }
        }
      }
    }
    my_argv = new char * [arg_vector.size()];
    for (size_t idx = 0; idx < arg_vector.size(); idx++)
    {
      assert(not arg_vector[idx].empty());
      if (replace_single_hyphen_options &&
          arg_vector[idx].size() > 2 &&
          arg_vector[idx][0] == '-' &&
          arg_vector[idx][1] != '-')
      {
        arg_vector[idx].insert(0,"-");
      }
      my_argv[idx] = strdup(arg_vector[idx].c_str());
    }
  }
  ~ArgvConcat(void)
  {
    for (size_t idx = 0; idx < arg_vector.size(); idx++)
    {
      free(my_argv[idx]);
    }
    delete[] my_argv;
  }
  void show() const noexcept
  {
    for (auto &&arg : arg_vector)
    {
      std::cout << "arg=" << arg << '\n';
    }
  }
  [[nodiscard]] size_t size() const noexcept { return arg_vector.size(); }
  [[nodiscard]] char **concat_get() const noexcept { return my_argv; }
};
#endif

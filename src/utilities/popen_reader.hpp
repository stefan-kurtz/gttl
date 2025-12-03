#ifndef POPEN_READER_HPP
#define POPEN_READER_HPP
#ifdef _WIN32
  #define NOMINMAX
  #include <io.h>
  #include "utilities/windows_fork.hpp"
#else
  #include <unistd.h>
  #include <stdio.h> // NOLINT(modernize-deprecated-headers)
#endif
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <cassert>

class PopenReader
{
  FILE *stdout_ptr,
       *stderr_ptr;
  public:
  /* allow to call the constructor with a variable number of arguments
     supplied to execlp. */
  template<typename... Args>
  PopenReader(const std::vector<std::string> &progs,Args ...args)
    : stdout_ptr(nullptr)
    , stderr_ptr(nullptr)
  {
    int stdout_tab[2];
    int stderr_tab[2];
#ifdef _WIN32
    int pipe_error = _pipe(stdout_tab, 4096UL * 1024 * 1024, 0);
#else
    int pipe_error = pipe(stdout_tab);
#endif
    if (pipe_error != 0)
    {
      fprintf(stderr,"failed to open stdout-pipe\n");
      exit(EXIT_FAILURE);
    }
#ifdef _WIN32
    pipe_error = _pipe(stderr_tab, 4096UL * 1024 * 1024, 0);
#else
    pipe_error = pipe(stderr_tab);
#endif
    if (pipe_error != 0)
    {
      fprintf(stderr,"failed to open stderr-pipe\n");
      exit(EXIT_FAILURE);
    }
    const pid_t pid = fork();
    if (pid == static_cast<pid_t>(0))
    {
#ifdef _WIN32
      dup2(stdout_tab[1], _fileno(stdout));
#else
      dup2(stdout_tab[1], STDOUT_FILENO);
#endif
      close(stdout_tab[0]);
      close(stdout_tab[1]);

#ifdef _WIN32
      dup2(stderr_tab[1], _fileno(stderr));
#else
      dup2(stderr_tab[1], STDERR_FILENO);
#endif
      close(stderr_tab[0]);
      close(stderr_tab[1]);

      for (auto const &prog : progs)
      {
        execlp(prog.c_str(), args..., (char *) nullptr);
      }
      exit(1);
    } else
    {
      if (pid < static_cast<pid_t>(0))
      {
        fprintf(stderr,"Child not created.\n");
        exit(EXIT_FAILURE);
      } else
      {
        close(stdout_tab[1]);
        close(stderr_tab[1]);
        stdout_ptr = fdopen(stdout_tab[0], "rb");
        stderr_ptr = fdopen(stderr_tab[0], "rb");
        return;
      }
    }
    exit(1);
  }
  size_t fread_stdout(void *ptr, size_t size, size_t count)
  {
    return std::fread(ptr, size, count, stdout_ptr);
  }
  int fgetc_stderr(void)
  {
    return std::fgetc(stderr_ptr);
  }
};
#endif

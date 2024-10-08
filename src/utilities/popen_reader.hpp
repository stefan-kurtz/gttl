#ifndef POPEN_READER_HPP
#define POPEN_READER_HPP
#ifdef _WIN32
  #define NOMINMAX
  #include <io.h>
#else
  #include <unistd.h>
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
    bool pipe_error = false;
    int stdout_tab[2], stderr_tab[2];
    pipe_error = static_cast<bool>(pipe(stdout_tab));
    pipe_error = pipe_error || static_cast<bool>(pipe(stderr_tab));
    assert(pipe_error == false);
    const pid_t pid = fork();
    if (pid == static_cast<pid_t>(0))
    {
      dup2(stdout_tab[1], STDOUT_FILENO);
      close(stdout_tab[0]);
      close(stdout_tab[1]);

      dup2(stderr_tab[1], STDERR_FILENO);
      close(stderr_tab[0]);
      close(stderr_tab[1]);

      for (auto const &prog : progs)
      {
        execlp(prog.c_str(), args..., (char *) NULL);
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
    return fread(ptr, size, count, stdout_ptr);
  }
  int fgetc_stderr(void)
  {
    return std::fgetc(stderr_ptr);
  }
};
#endif

#ifndef WINDOWS_MKDTEMP_HPP
#define WINDOWS_MKDTEMP_HPP
#ifdef _WIN32
#include <cstdio>
#include <cstring>
#include <direct.h>
#define NOMINMAX
#include <io.h>
#include <cstdlib>

char* mkdtemp(char* template_str)
{
  if(_mktemp_s(template_str, strlen(template_str) + 1) != 0)
  {
    return nullptr;
  }
  if(_mkdir(template_str) != 0)
  {
    return nullptr;
  }
  return template_str;
}

int rmdir(const char* dirname)
{
  return _rmdir(dirname);
}

#endif
#endif // WINDOWS_MKDTEMP_HPP

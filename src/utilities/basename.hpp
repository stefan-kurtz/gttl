#ifndef BASENAME_HPP
#define BASENAME_HPP
#include <cstring>
#include <cstdlib>

class GttlBasename
{
 private:
  char *sbuf;

 public:
  GttlBasename(const char *path)
  {
    char *c;
    bool foundother = false;
    size_t i, pathlen;

    if (path != NULL)
    {
      pathlen = strlen(path);
    } else
    {
      pathlen = 0;
    }
    sbuf = static_cast<char *>(malloc((pathlen + 2) * sizeof(char)));
    if (path == NULL || *path == '\0')
    {
      strcpy(sbuf, ".");
      return;
    }
    strcpy(sbuf, path);
    for (c = sbuf + pathlen - 1; c >= sbuf; c--)
    {
      if (*c == '/')
      {
        if (foundother)
        {
          c++;
          for (i = 0; c[i] != '\0'; i++)
          {
            sbuf[i] = c[i];
          }
          sbuf[i] = '\0';
          break;
        }
        if (c > sbuf)
        {
          *c = '\0';
        }
      } else
      {
        foundother = true;
      }
    }
  }
  const char *str(void) { return sbuf; }
  ~GttlBasename(void) { free(sbuf); }
};
#endif

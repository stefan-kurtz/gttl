# Portability

in order to maintain portability between operating systems, the following rules
should be followed:

## Printf-Specifiers

All printf/fprintf/etc. calls should use the following format-specifiers to
print the given  data-types:

| Datatype            | Format-specifier |
|:-------------------:|:----------------:|
| size_t              | %zu              |
| unsigned int        | %u               |
| unsigned long (int) | %lu              |
| unsigned long long  | %llu             |
| int                 | %d               |
| long (int)          | %ld              |
| long long           | %lld             |
| unsigned char       | %u               |
| signed char         | %d               |


## Makefiles and Scripts

- All calls to `diff` need to include the `--strip-trailing-cr` flag
- All calls to `mktemp` need to include the `--tmpdir=.` argument, since
`/tmp` is not a directory on Windows.
    - The `make clean`-target thus needs to remove the created temporary files.
- `egrep` and `fgrep` may not exist, use `grep -E` and `grep -F` respectively.
- `CFLAGS`, `LDLIBS`, etc. should be copied according to the patterns used
in `testsuite/Makefile`.

## Includes

We can tell whether we are on a Windows-machine at compile-time by using the
`#ifdef _WIN32` preprocessor-macro.

- `<unistd.h>` does not exist on Windows. Most of its features can be found in
`<io.h>`.
    - When including `<io.h>`, we must first `#define NOMINMAX` to not break `std::min`/`std::max`.
- `<sys/mman.h>` does not exist on Windows. The relevant features used in GTTL
can be found in `"utilities/windows_mman.hpp"`.
- In order to use `getopt` and its related variables,
`"utilities/windows_getopt.hpp"` must be included on Windows.

## Dependencies

ZLib is an external dependency that needs to be explicitly installed on
Windows, eg. with
```
vcpkg install zlib
```

## Gitattributes

All files should have their line-end conversions (or lack thereof) explicitly
defined in a `.gitattributes` file. This is verified by `make code_check`.

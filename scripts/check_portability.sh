#!/bin/bash

errors_exist=false

git --no-pager grep -e "diff " --and \( --not -e "--strip-trailing-cr" \) -- '**Makefile' '**.sh' '**.py' && \
    printf "%s: diff should only be used with the --strip-trailing-cr flag!\n\n" "$0" && \
    errors_exist=true

git --no-pager grep -e "%lu" -- '**.cpp' '**.hpp' && \
    printf "%s: size_t variables should be used where possible and should be printed using %%zu.\n\n" "$0" && \
    errors_exist=true

git --no-pager grep -e "mktemp " --and \( --not -e "--tmpdir=" \) -- '**Makefile' '**.sh' '**.py' && \
    printf "%s: mktemp needs to be called with --tmpdir=. or similar, since /tmp does not exist on Windows.\n\n" "$0" && \
    errors_exist=true

git --no-pager grep -e "egrep" -e "fgrep" -- '**Makefile' '**.sh' '**.py' ':(exclude)scripts/check_portability.sh' && \
    printf "%s: egrep and fgrep may not exist on all systems. Use grep -E and grep -F instead.\n\n" "$0" && \
    errors_exist=true

git --no-pager grep -L -e "#define NOMINMAX" $(git grep -l -e "#include <io.h>") && \
    printf "%s: You should #define NOMINMAX before including <io.h> in the above files!\n\n" "$0" && \
    errors_exist=true

if [ "$errors_exist" = true ] ; then
    exit 1
fi

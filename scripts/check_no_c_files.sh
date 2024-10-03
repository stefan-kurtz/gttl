#!/bin/bash

# Checks that no .c or .h files are found within the repository.
# Exits with return-code 1 if any such files are found, 0 otherwise.
# Used as part of `make code_check`.

C_FILES=$(find . -type f \( -name '*.c' -o -name '*.h' \))

if [ -n "${C_FILES}" ]; then
  >&2 echo "The following C-Files were found in the repository:"
  for c_file in "${C_FILES}"; do
    >&2 echo "${c_file}"
  done
  exit 1
fi

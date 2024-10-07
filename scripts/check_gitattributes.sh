#!/bin/bash

# Checks that all files in the current git-repository are mentioned in
# .gitattributes by checking whether their `text` attribute was set or unset
# Exits with return-code 1 if any checked-in files are not missing, 0 otherwise.
# This automatically ignores all files in `.gitignore`.
# Used as part of `make code_check`.

FILES_WITHOUT_ATTR=$(git ls-files | xargs git check-attr text |\
  egrep -v "text: set|text: unset")

if [ -n "${FILES_WITHOUT_ATTR}" ]; then
  >&2 echo "The following Files were found in the repository, but are not mentioned in .gitattributes:"
  for file in "${FILES_WITHOUT_ATTR}"; do
    >&2 echo "${file}"
  done
  exit 1
fi

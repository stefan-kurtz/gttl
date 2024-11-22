#!/bin/sh
valgrind --quiet --tool=memcheck --memcheck:leak-check=full --memcheck:leak-resolution=high --error-exitcode=1 --log-fd=1 --error-limit=yes --dsymutil=yes ./"$@"

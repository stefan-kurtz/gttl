#!/bin/sh

if test $# -ne 1
then
  echo "Usage: $0 <inputfile>"
  exit 1
fi

set -e -x

inputfile=$1
./sa_induced.x --verbose --indexname sa --lcptab plcp5n --absolute_suftab ${inputfile}
./sa_induced.x --verbose --indexname sa --succinct --lcptab plcp5n --absolute_suftab ${inputfile}
./lcp_checker.x sa
for suffix in prj lcp ll2 ll4 lls suf tis isa bsf
do
  rm -f sa.${suffix}
done

#!/bin/sh

set -e -x

if test $# -ne 1
then
  echo "Usage: $0 <inputfile>"
  exit 1
fi

inputfile=$1
TMPFILE=`mktemp TMP.XXXXXX` || exit 1
./char_range.py --alphabet 'N' ${inputfile} > ${TMPFILE}
./char_range_mn.x --singlechar ${inputfile} | diff --strip-trailing-cr - ${TMPFILE}
for opti in '' --invert
do
  for optr in '' --reverse
  do
    for optm in '' --multiseq
    do
      ./char_range.py ${opti} ${optr} ${optm}  ${inputfile} > ${TMPFILE}
      ./char_range_mn.x ${opti} ${optr} ${optm} ${inputfile} | diff --strip-trailing-cr - ${TMPFILE}
    done
  done
done
rm -f ${TMPFILE}

#!/bin/sh

if test $# -ne 1
then
  echo "Usage: $0 <inputfile>"
  exit 1
fi

inputfile=$1
TMPFILE1=`mktemp --tmpdir=.` || exit 1
./multiseq_mn.x --width 70 --sorted_by_header ${inputfile} > ${TMPFILE1}
TMPFILE2=`mktemp --tmpdir=.` || exit 1
grep -v '^#' ${TMPFILE1} | grep '^>' > ${TMPFILE2}
grep '^>' ${inputfile} | sort | diff --strip-trailing-cr - ${TMPFILE2}
../scripts/cmp_sorted_sequences.py ${inputfile} ${TMPFILE1}
rm -f ${TMPFILE1} ${TMPFILE2}

#!/bin/sh

if test $# -ne 1
then
  echo "Usage: $0 <inputfile>"
  exit 1
fi

inputfile=$1
TMPFILE1=`mktemp TMP.XXXXXX` || exit 1
./multiseq_mn.x --width 70 --sorted_by_header ${inputfile} > ${TMPFILE1}
TMPFILE2=`mktemp TMP.XXXXXX` || exit 1
grep -v '^#' ${TMPFILE1} | grep '^>' > ${TMPFILE2}
grep '^>' ${inputfile} | sort | diff - ${TMPFILE2}
../scripts/cmp_sorted_by_header.py ${inputfile} ${TMPFILE1}
rm -f ${TMPFILE1} ${TMPFILE2}

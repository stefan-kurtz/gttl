#!/bin/sh

set -e -x

TMPFILE1=`mktemp --tmpdir=. TMP.XXXXXX` || exit 1
TMPFILE2=`mktemp --tmpdir=. TMP.XXXXXX` || exit 1
TMPFILE3=`mktemp --tmpdir=. TMP.XXXXXX` || exit 1
TMPFILE4=`mktemp --tmpdir=. TMP.XXXXXX` || exit 1
echo "atggc" > ${TMPFILE1}
echo ">\naggtatcgc" > ${TMPFILE2}
echo "atggc	5" > ${TMPFILE3}
for strategy in p_times_s s_times_p
do
  ./myersapm.x ${strategy} 1 ${TMPFILE1} ${TMPFILE2} 2 > ${TMPFILE4}
  diff -I '^#' --strip-trailing-cr ${TMPFILE3} ${TMPFILE4}
done
rm -f ${TMPFILE1} ${TMPFILE2} ${TMPFILE3} ${TMPFILE4}

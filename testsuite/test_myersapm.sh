#!/bin/sh

set -e -x

TMPFILE1=`mktemp --tmpdir=. TMP.XXXXXX` || exit 1
TMPFILE2=`mktemp --tmpdir=. TMP.XXXXXX` || exit 1
TMPFILE3=`mktemp --tmpdir=. TMP.XXXXXX` || exit 1
TMPFILE4=`mktemp --tmpdir=. TMP.XXXXXX` || exit 1
echo "atggc" > ${TMPFILE1}
echo ">\naggtatcgc" > ${TMPFILE2}
./myersapm.x ${TMPFILE1} ${TMPFILE2} 2 > ${TMPFILE3}
echo "atggc	5" > ${TMPFILE4}
diff ${TMPFILE3} ${TMPFILE4}
rm -f ${TMPFILE1} ${TMPFILE2} ${TMPFILE3} ${TMPFILE4}

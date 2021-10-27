#!/bin/sh

set -e -x

if test $# -ne 1
then
  echo "Usage: $0 <inputfile>"
  exit 1
fi

inputfile=$1
TMPFILE=`mktemp TMP.XXXXXX` || exit 1
./char_range.py ${inputfile} > ${TMPFILE}
./char_range_mn.x ${inputfile} | diff - ${TMPFILE}
rm -f ${TMPFILE}

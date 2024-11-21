#!/bin/sh

if test $# -ne 3
then
  echo  "Usage: $0 <cutlen> <rows> <cols>"
  exit 1
fi

cutlen=$1
rows=$2
cols=$3
TMPFILE=`mktemp --tmpdir=.` || exit 1
./matrix_partition.x ${cutlen} ${rows} ${cols} > ${TMPFILE}
./matrix_partition.py ${cutlen} ${rows} ${cols} | diff --strip-trailing-cr - ${TMPFILE}
rm -f ${TMPFILE}

#!/bin/sh

if test $# -ne 1
then
  echo "Usage: $0 <inputfile>"
  exit 1
fi

inputfile=$1
TMPFILE=`mktemp --tmpdir=.` || exit 1
./multiseq_mn.x --width 70 --sorted_by_length ${inputfile} > ${TMPFILE}
../scripts/cmp_sorted_sequences.py --sorted_by_length ${inputfile} ${TMPFILE}
rm -f ${TMPFILE}

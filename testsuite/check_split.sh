#!/bin/sh

set -e

if test $# -ne 1
then
  echo "Usage: $0 <inputfile in fastq_format"
  exit 1
fi

inputfile=$1
inputfile_basename=`basename ${inputfile}`
for split_size in 13 29 42 79
do
  ./fastq_mn.x -e --split_size ${split_size} ${inputfile}
  cat `ls ${inputfile_basename}*` | diff - ${inputfile}
  rm -f `ls ${inputfile_basename}*`
done
./fastq_mn.x -e ${inputfile} | diff - ${inputfile}

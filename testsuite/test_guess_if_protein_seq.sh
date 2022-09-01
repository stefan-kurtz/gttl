#!/bin/sh 

for filename in `ls ../testdata/*.fna ../testdata/*.fasta`
do
  grep -q `basename ${filename}` ../testdata/corrupt_files.txt
  if test $? -eq 1
  then
    ../scripts/guess_if_protein_seq.py ${filename}
    result1=$?
    ./guess_if_protein_seq.x ${filename}
    result2=$?
    if test ${result1} -ne ${result2}
    then
      echo "$0: inconsistent guesses for ${filename}"
      exit 1
    else
      if test ${result1} -eq 0
      then
        echo "${filename} contains protein sequences"
      else
        echo "${filename} contains DNA sequences"
      fi
    fi
  fi
done

#!/bin/sh

set -e -x

outname=test_splitseq
cat_file=tmp_complete.fasta
file_list=`ls ../../testdata/vaccg.fna ../../testdata/at1MB.fna ../../testdata/ychrIII.fna`

cat ${file_list} > ${cat_file}

options="--lcptab kasai13n --check_suftab --relative_suftab"

./sa_induced.x ${options} --indexname ${outname} ${cat_file}
./sa_induced.x ${options} --indexname ${outname}_parts ${file_list}

for suffix in bsf lcp ll2 ll4
do
  diff --strip-trailing-cr ${outname}.${suffix} ${outname}_parts.${suffix}
done
rm -f ${outname}.* ${outname}_parts.* ${cat_file}

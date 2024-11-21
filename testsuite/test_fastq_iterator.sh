#!/bin/sh

set -e -x

for fastqfile in 70x_161nt_phred64.fastq SRR19536726_1_1000.fastq.gz varlen_paired_2.fastq varlen_paired_2.fastq
do
  for bits in 8 16 32 64
  do
    ./fastq_mn.x --encoding ${bits} ../testdata/${fastqfile}
  done
done

for fastqfile in 70x_161nt_phred64.fastq SRR19536726_1_1000.fastq.gz
do
  qgram_length=2
  while test ${qgram_length} -lt 32
  do
    ./fastq_mn.x --encoding 64,${qgram_length} ../testdata/${fastqfile}
    qgram_length=`expr ${qgram_length} + 1`
  done
done

for mapped_opt in "" --mapped
do
  TMPFILE=`mktemp --tmpdir=.` || exit 1
  touch ${TMPFILE}
  for fastqfile in `ls ../testdata/*.fastq | sort`
  do
    fastqfile=$(echo ${fastqfile} | tr -d '\r')
    ./fastq_mn.x ${mapped_opt} --echo $fastqfile | diff --strip-trailing-cr - $fastqfile
    ./fastq_mn.x ${mapped_opt} --hash wy ${fastqfile} >> ${TMPFILE}
  done
  diff --strip-trailing-cr ${TMPFILE} ../testdata/wy_hash_values.tsv
  rm -f ${TMPFILE}
  ./fastq_mn.x ${mapped_opt} --statistics ../testdata/varlen_paired_1.fastq | \
     sed -e 's/paired_1/paired/' | diff --strip-trailing-cr - ../testdata/varlen_paired_stat.tsv
  ./fastq_mn.x ${mapped_opt} --statistics ../testdata/varlen_paired_2.fastq | \
     sed -e 's/paired_2/paired/' | diff --strip-trailing-cr - ../testdata/varlen_paired_stat.tsv
  ./fastq_mn.x ${mapped_opt} --statistics ../testdata/70x_161nt_phred64.fastq |\
      diff --strip-trailing-cr - ../testdata/70x_161nt_phred64_stat.tsv
  ./fastq_mn.x ${mapped_opt} --fasta_output ../testdata/varlen_paired_1.fastq |\
      diff --strip-trailing-cr - ../testdata/varlen_paired_1.fasta
done
./fastq_mn.x --fasta_output ../testdata/varlen_paired_[12].fastq | diff --strip-trailing-cr - ../testdata/varlen_paired_both.fasta

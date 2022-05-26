#!/bin/sh

set -e -x

for mapped_opt in "" --mapped 
do
  for fastqfile in `ls ../testdata/*.fastq`
  do
   ./fastq_mn.x ${mapped_opt} --echo $fastqfile | diff - $fastqfile
  done
  ./fastq_mn.x ${mapped_opt} --statistics ../testdata/varlen_paired_1.fastq | diff - ../testdata/varlen_paired_stat.tsv
  ./fastq_mn.x ${mapped_opt} --statistics ../testdata/varlen_paired_2.fastq | diff - ../testdata/varlen_paired_stat.tsv
  ./fastq_mn.x ${mapped_opt} --statistics ../testdata/70x_161nt_phred64.fastq | diff - ../testdata/70x_161nt_phred64_stat.tsv
  ./fastq_mn.x ${mapped_opt} --fasta_output ../testdata/varlen_paired_1.fastq | diff -  ../testdata/varlen_paired_1.fasta
done
./fastq_mn.x --fasta_output ../testdata/varlen_paired_[12].fastq | diff - ../testdata/varlen_paired_both.fasta

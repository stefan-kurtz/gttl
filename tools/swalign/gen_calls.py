#!/usr/bin/env python3

import sys, os
from wheels import cartproductsingleloop

def combiner(lists):
  set_size_tab = list(map(len,lists))
  for wheel in cartproductsingleloop(set_size_tab):
    args = ' '.join([lists[idx][choice] for idx, choice in enumerate(wheel)])
    print('./sw_all_against_all.x {}'.format(args))

input_files = ['seq18.fasta',
               'smrt100-2000-first.fna',
               'query-protein_f5.fsa',
               'smrt100-2000-last.fna']

for filename in input_files:
  filepath = 'testdata/{}'.format(filename)
  if not os.path.isfile(filepath):
    sys.stderr.write('{}: file {} does not exist\n'.format(filepath))
    exit(1)

print('#!/bin/sh')
print('set -e -x')

input_cmds = ['-d testdata/seq18.fasta -c 50',
              '-d testdata/seq18.fasta',
              '-d testdata/seq18.fasta -s unit_score_aa',
              '-d testdata/smrt100-2000-first.fna',
              '-d testdata/seq18.fasta -q testdata/query-protein_f5.fsa -c 50',
              '-d testdata/seq18.fasta -q testdata/query-protein_f5.fsa',
              '-d testdata/seq18.fasta -q testdata/query-protein_f5.fsa '
                  '-s unit_score_aa',
              '-d testdata/smrt100-2000-first.fna -q '
                  'testdata/smrt100-2000-last.fna',
              '-d testdata/query-protein_f5.fsa -q testdata/seq18.fasta',
              '-d testdata/smrt100-2000-last.fna '
                  '-q testdata/smrt100-2000-first.fna']

gap_scoring = ['-g 7 2','']
vectorized = ['-v {}'.format(v) for v in [1,2]]
output = ['-a {}'.format(a) for a in [1,2,60]]
threads = ['-t {}'.format(t) for t in [1,2,4]]
minimize_space = ['-m','']
header_display = ['-h','']
best_mode = ['-b 2','']

combiner([input_cmds,gap_scoring,vectorized,output,\
          threads,minimize_space,header_display,best_mode])

#!/usr/bin/env python3

import sys, argparse
from fastaIterator import fasta_next

def parse_arguments(argv):
  p = argparse.ArgumentParser(description=('enumerate maximum ranges of '
                                           'characters from some given '
                                           'alphaber'))
  p.add_argument('inputfile',type=str,
                  help='specify input file')
  return p.parse_args(argv)

def enum_char_ranges(seq,alphabet):
  previous_start = None
  previous_length = None
  for idx, cc in enumerate(seq):
    #print('idx={},cc={}'.format(idx,cc))
    if previous_length is None:
      if cc in alphabet:
        previous_start = idx
        previous_length = 1
    else:
      assert previous_start is not None
      if cc in alphabet:
        previous_length += 1
      else:
        yield previous_start, previous_length
        previous_length = None
        previous_start = None
  if previous_length is not None:
    yield previous_start,previous_length

args = parse_arguments(sys.argv[1:])

non_wildcard_ranges_total_length = 0
for seqnum, (header, sequence) in enumerate(fasta_next(args.inputfile)):
  for start, length in enum_char_ranges(sequence,set('ACGTacgt')):
    print('{}\t{}\t{}'.format(seqnum,start,length))
    non_wildcard_ranges_total_length += length
print('# non_wildcard_ranges_total_length\t{}'
       .format(non_wildcard_ranges_total_length))

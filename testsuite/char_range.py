#!/usr/bin/env python3

import sys, argparse
from fastaIterator import fasta_next

def parse_arguments(argv):
  p = argparse.ArgumentParser(description=('enumerate maximum ranges of '
                                           'characters which are/are not '
                                           'from some given alphabet'))
  p.add_argument('-i','--invert',action='store_true',default=False,
                  help='inversion, i.e. are not from alphabet')
  p.add_argument('-r','--reverse',action='store_true',default=False,
                  help='enumerate ranges in reverse order')
  p.add_argument('inputfiles',nargs='+',help='specify input files')
  return p.parse_args(argv)

def enum_char_ranges(invert,seq,alphabet):
  previous_start = None
  previous_length = None
  for idx, cc in enumerate(seq):
    #print('idx={},cc={}'.format(idx,cc))
    condition = cc in alphabet
    if invert:
      condition = not condition
    if previous_length is None:
      if condition:
        previous_start = idx
        previous_length = 1
    else:
      assert previous_start is not None
      if condition:
        previous_length += 1
      else:
        yield previous_start, previous_length
        previous_length = None
        previous_start = None
  if previous_length is not None:
    yield previous_start,previous_length

args = parse_arguments(sys.argv[1:])

for inputfile in args.inputfiles:
  print(inputfile)
  non_wildcard_ranges_total_length = 0
  for seqnum, (header, sequence) in enumerate(fasta_next(inputfile)):
    r_list = list()
    for start, length in enum_char_ranges(args.invert,sequence,set('ACGTacgt')):
      non_wildcard_ranges_total_length += length
      if args.reverse:
        r_list.append((start,length))
      else:
        print('{}\t{}\t{}'.format(seqnum,start,length))
    if args.reverse:
      for start, length in reversed(r_list):
        print('{}\t{}\t{}'.format(seqnum,start,length))
  print('# non_wildcard_ranges_total_length\t{}'
         .format(non_wildcard_ranges_total_length))

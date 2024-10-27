#!/usr/bin/env python3

import sys, argparse
from fastaIterator import fasta_next

def parse_command_line(argv):
  p = argparse.ArgumentParser(description=('verify that sequences are '
                                           'correctly sorted by headers and '
                                           'that no sequence information is '
                                           'lost'))
  p.add_argument('origfile',type=str,help='specify path of original file')
  p.add_argument('sortedfile',type=str,help='specify path of sorted file')
  return p.parse_args(argv)

def collect_sequences(inputfile):
  seq_dict = dict()
  for header, sequence in fasta_next(inputfile):
    if header in seq_dict:
      sys.stderr.write('{}: duplicate header {}\n'.format(sys.argv[0],header))
      exit(1)
    seq_dict[header] = sequence
  return seq_dict

args = parse_command_line(sys.argv[1:])

orig_dict = collect_sequences(args.origfile)
sorted_headers = list(sorted(orig_dict))
checked = 0
for idx, (header, sequence) in enumerate(fasta_next(args.sortedfile)):
  if header != sorted_headers[idx]:
    sys.stderr.write('{}: header {} != {}\n'.format(sys.argv[0],header,
                                                    sorted_headers))
    exit(1)
  if header not in orig_dict:
    sys.stderr.write('{}: header {} does not occur in sequences of file {}\n'
                     .format(sys.argv[0],header,
                                                      orig_dict[header]))
    exit(1)
  if sequence != orig_dict[header]:
    sys.stderr.write('{}: header {} with sequence(length{}): {} != {}\n'
                     .format(sys.argv[0],header,len(sequence),
                             sequence,orig_dict[header]))
    exit(1)
  checked += 1
assert checked == len(orig_dict)

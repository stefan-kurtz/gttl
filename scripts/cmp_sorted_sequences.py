#!/usr/bin/env python3

import sys, argparse
from fastaIterator import fasta_next

def parse_command_line(argv):
  p = argparse.ArgumentParser(description=('verify that sequences are '
                                           'correctly sorted and that no '
                                           'sequence information is lost'))
  p.add_argument('-l','--sorted_by_length',action='store_true',default=False,
                 help='sort the sequences by length; default: sort by header')
  p.add_argument('origfile',type=str,help='specify path of original file')
  p.add_argument('sortedfile',type=str,help='specify path of sorted file')
  return p.parse_args(argv)

def collect_sequences(inputfile):
  seq_dict = dict()
  for header, sequence in fasta_next(inputfile):
    if header in seq_dict:
      sys.stderr.write(f'{sys.argv[0]}: duplicate header {header}\n')
      exit(1)
    seq_dict[header] = sequence
  return seq_dict

args = parse_command_line(sys.argv[1:])

orig_dict = collect_sequences(args.origfile)
if args.sorted_by_length:
  sorted_headers = sorted(orig_dict, key=lambda h: len(orig_dict[h]))
else:
  sorted_headers = sorted(orig_dict)
checked = 0
for idx, (header, sequence) in enumerate(fasta_next(args.sortedfile)):
  if header != sorted_headers[idx]:
    sys.stderr.write(f'{sys.argv[0]}: header {header} != {sorted_headers}\n')
    exit(1)
  if header not in orig_dict:
    sys.stderr.write(f'{sys.argv[0]}: header {header} does not occur '
                     f'in sequences of file {orig_dict[header]}\n')
    exit(1)
  if sequence != orig_dict[header]:
    sys.stderr.write(f'{sys.argv[0]}: header {header} with '
                     f'sequence(length {len(sequence)}): '
                     f'{sequence} != {orig_dict[header]}\n')
    exit(1)
  checked += 1
assert checked == len(orig_dict)

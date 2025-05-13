#!/usr/bin/env python3

import sys, argparse
from fastaIterator import fasta_next

def parse_arguments(argv):
  p = argparse.ArgumentParser(description=('decide if input file contains '
                                           'protein sequences'))
  p.add_argument('-p','--print_result',action='store_true',default=False,
                  help=('show result by printing 1 (for true) or 0 (for false) '
                        ' and always exit with code 0, if no other problems '
                        'occur'))
  p.add_argument('inputfiles',nargs='+',
                  help='specify input files, - means stdin')
  return p.parse_args(argv)

def guess_if_protein_sequence(guess_size_to_decide,sequence):
  protein_only_characters = set(list('LIFEQPXZ'))
  for _, cc in zip(range(guess_size_to_decide),sequence):
    if cc in protein_only_characters:
      return True
  return False

def guess_if_protein_file(inputfiles):
  total_length = 0
  decided_if_protein = False
  guess_size_to_decide = 1000
  stop = False
  for inputfile in inputfiles:
    for _, sequence in fasta_next(inputfile):
      if guess_if_protein_sequence(guess_size_to_decide,sequence):
        decided_if_protein = True
        stop = True
        break
      total_length += len(sequence)
      if total_length >= guess_size_to_decide:
        stop = True
        break
    if stop:
      break
  return decided_if_protein

args = parse_arguments(sys.argv[1:])

is_protein = guess_if_protein_file(args.inputfiles)
if args.print_result:
  print(1 if is_protein else 0)
else:
  if is_protein:
    exit(0)
  exit(1)

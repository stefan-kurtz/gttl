#!/usr/bin/env python3
import sys, argparse
from fastaIterator import fasta_next

def parse_command_line(argv):
  p = argparse.ArgumentParser(description='check correctness of unwords')
  p.add_argument('unwords_file',type=str,
                  help=('specify file with unwords, one per line not beginning '
                        'with #'))
  p.add_argument('sequence_file',type=str,
                  help=('specify file with sequences from which the unwords '
                        'were computed'))
  return p.parse_args(argv)


args = parse_command_line(sys.argv[1:])

# Vorverarbeitung der unwords
try:
  unwords_stream = open(args.unwords_file)
except IOError as err:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
  exit(1)

unwords_list = [unword.rstrip() for unword in unwords_stream \
                                if not unword.startswith('#')]
unwords_stream.close()

# sequenzweise "Uberpr"ufung der unwords
no_unwords_found = True
for _,seq in fasta_next(args.sequence_file):
  seq = seq.upper()
  for unword in unwords_list:
    if seq.find(unword) != -1:
      no_unwords_found = False
      print('{} is incorrect.'.format(unword))

if no_unwords_found:
  print('All unwords are correct.')

#!/usr/bin/env python3

import os

def enum_big_fasta_files(environment_vars):
  for evar in environment_vars:
    if evar in os.environ:
      filepath = os.environ[evar]
      if os.path.isfile(filepath):
        yield filepath

if __name__ == '__main__':
  import sys, argparse
  p = argparse.ArgumentParser(description=('output list of files files from a '
                                           'given list of environment vars'))
  p.add_argument('environ_vars',nargs='+',
                  help='specify enivorment vars')
  args = p.parse_args(sys.argv[1:])
  print('\n'.join(enum_big_fasta_files(args.environ_vars)))

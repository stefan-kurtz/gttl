#!/usr/bin/env python3

import os

def enum_gttl_fasta_files(directory):
  file_list_fasta = ['at1MB.fna','simple.fna','small.fna','vaccg.fna',
                     'wc_palindromes8.fna','ychrIII.fna','sw175.fna']
  path_list_fasta = [f'{directory}/{filename}' for filename in file_list_fasta]
  for p in path_list_fasta:
    assert os.path.isfile(p), f'{sys.argv[0]}: {p} does not exist\n'
  return path_list_fasta

if __name__ == '__main__':
  import sys, argparse
  p = argparse.ArgumentParser(description=('output list of fasta files from '
                                           'gttl/testdata'))
  p.add_argument('directory',type=str,
                  help='specify input file, means stdin')
  args = p.parse_args(sys.argv[1:])
  print('\n'.join(enum_gttl_fasta_files(args.directory)))

#!/usr/bin/env python3

import re, os, sys, argparse, subprocess, shlex, pathlib
sys.path.append('../../scripts')
from guess_if_protein_seq import guess_if_protein_file

def parse_arguments(argv):
  p = argparse.ArgumentParser(description='run sa_induced.x test cases')
  p.add_argument('-d','--dry_run',action='store_true',default=False,
                  help='show commands instead of executing them')
  p.add_argument('-a','--absolute_suftab',action='store_true',default=False,
                  help='output absolute suftab')
  p.add_argument('-r','--relative_suftab',action='store_true',default=False,
                  help='output relative suftab')
  p.add_argument('-l','--lcptab',type=str,default=None,
                  choices=['kasai13n','kasai9n','plcp5n'],
                  help=('construct lcptable using one of the three possible '
                        'implementations, default: None'))
  p.add_argument('--succinct',action='store_true',default=False,
                 help=('also construct the succinct representation of '
                       'lcp-table in file with suffix .lls'))
  p.add_argument('--reverse_complement',action='store_true',default=False,
                 help='consider the reverse complement if sequence is DNA')
  p.add_argument('-v','--verbose',action='store_true',default=False,
                  help='run in verbose mode')
  return p.parse_args(argv)

def split_lines(bstring):
  return [bs.decode() for bs in bstring.split(b'\n')]

def is_index_file(filename, index_suffixes):
  for suffix in index_suffixes:
    if filename.endswith(suffix):
      return True
  return False

def my_run(cmd):
  try:
    subprocess.run(shlex.split(cmd))
  except subprocess.CalledProcessError as exc:
    sys.stderr.write(f'{sys.argv[0]}: {cmd} failed with returncode '
                     f'{exc.returncode} failed, output: {exc.output}\n')
    exit(1)

args = parse_arguments(sys.argv[1:])

index_suffixes = ['prj','lcp','ll2','ll4','lls','suf','tis','isa','bsf']

verbose_opt = '-v' if args.verbose else ''
file_list_plain_paths = [filename for filename in os.listdir('.')
                                  if os.path.getsize(filename) > 0
                                  if os.path.isfile(filename)
                                  if not is_index_file(filename,index_suffixes)]

if 'FASTA_FILES' in os.environ:
  fasta_files = os.environ.get('FASTA_FILES').split(' ')
else:
  fasta_files = [f'{os.environ["GTTL"]}/testdata/at1MB.fna']

if args.dry_run:
  print('#!/bin/sh')
  print('set -e -x')
for file_list, plain in [(file_list_plain_paths,True), (fasta_files,False)]:
  for filename in file_list:
    if filename == '':
      continue
    if args.reverse_complement and not plain and \
       not guess_if_protein_file([filename]):
      reverse_complement = '--reverse_complement'
    else:
      reverse_complement = ''
    if args.absolute_suftab or args.lcptab in ['kasai9n','plcp5n']:
      absolute_suftab = '--absolute_suftab'
    else:
      absolute_suftab = ''
    if args.lcptab is None or args.lcptab == 'kasai13n':
      check_suftab = '--check_suftab'
    else:
      check_suftab = ''
    option_list = [verbose_opt,
                   '--plain_input_format' \
                     if plain else '',
                   absolute_suftab,
                   '--relative_suftab' if (not plain and \
                                           args.relative_suftab) else '',
                   f'--lcptab {args.lcptab}' if args.lcptab else "",
                   check_suftab,
                   reverse_complement,
                   filename]
    clean_option_list = [opt for opt in option_list if opt != '']
    cmd = f'./sa_induced.x {" ".join(clean_option_list)}'
    indexname = os.path.basename(filename)
    if args.dry_run:
      print(cmd)
    else:
      my_run(cmd)
    if not plain and args.lcptab and args.succinct:
      cmd1 = f'./sa_induced.x --succinct {" ".join(clean_option_list)}'
      cmd2 = f'./lcp_checker.x {indexname}'
      for cmd in cmd1,cmd2:
        if args.dry_run:
          print(cmd)
        else:
          my_run(cmd)
  for suffix in index_suffixes:
    this_filename = f'{indexname}.{suffix}'
    if os.path.isfile(this_filename):
      pathlib.Path(this_filename).unlink()

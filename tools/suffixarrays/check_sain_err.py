#!/usr/bin/env python3

import argparse, sys, shutil, os
sys.path.append('../../scripts')
from mysubprocess import mysubprocess_expect

def parse_arguments():
  p = argparse.ArgumentParser(description='run program for error cases')
  p.add_argument('-e','--echo',type=str,default=None,
                  help='echo program calls with given prefix')
  p.add_argument('-v','--valgrind',action='store_true',default=False,
                  help='run program with valgrind')
  p.add_argument('path',type=str,default=None,
                  help='specify path of test files')
  return p.parse_args()

args = parse_arguments()

if args.valgrind and os.uname().sysname != 'Darwin' and \
                     shutil.which('valgrind'):
  valgrind_opts = ['--quiet', '--tool=memcheck', '--memcheck:leak-check=full',
                   '--memcheck:leak-resolution=high','--error-exitcode=1',
                   '--log-fd=1','--error-limit=yes','--dsymutil=yes']
  prefix = 'valgrind {} '.format(' '.join(valgrind_opts))
else:
  prefix = ''

'''
each call to the following dictionary entries represent a corrupted
file and the expected error message with placeholders for the
prefix and the filename, when applying the multiseq_test.x to that file.
the expected error code (which is EXIT_FAILURE in the
C/C++-source) is always 1.
'''

prot_file = '{}/testdata/sw175.fna'.format(args.path)
test_cases = [(('./sa_induced.x {}/testdata/non_existing.fna')
                .format(args.path),
               ('./sa_induced.x: file "{}/testdata/non_existing.fna": '
                'cannot open file')
                .format(args.path,args.path)),
               ('./sa_induced.x --plain_input_format {}/testdata/empty.fna'
                .format(args.path),
                ('./sa_induced.x: file "{}/testdata/empty.fna"'
                 ': cannot construct suffix array from empty string')
                 .format(args.path)),
                ('./sa_induced.x --reverse_complement {}'.format(prot_file),
                 ('./sa_induced.x: file "{}": option --reverse_complement is '
                  'only possible for DNA sequences').format(prot_file))
              ]

dna_file = '{}/testdata/at1MB.fna'.format(args.path)
for plaininput in [True,False]:
  for option, suffix in [('-a','suf'),('-r','bsf')]:
    test_cases.append(('./sa_induced.x -o unwritable_dir/tmp {} {} {}'\
                       .format(option,
                               '--plain_input_format' \
                                  if plaininput and suffix == 'suf' else '',
                               dna_file),
                       ('./sa_induced.x: file "{}": cannot create file '
                       '"unwritable_dir/tmp.{}"').format(dna_file,suffix)))

if args.echo:
  print('#!/bin/sh')
  print('set -e -x')
expected_error_code = 1
for call, expected_err_msg in test_cases:
  prefix_call = '{}{}'.format(prefix,call)
  if args.echo:
    print('{} {}'.format(args.echo,prefix_call))
  else:
    mysubprocess_expect(prefix_call,expected_error_code,expected_err_msg)

if not args.echo:
  sys.stderr.write('Congratulations: {} checked {} test cases for ./sa_induced.x\n'
                    .format(sys.argv[0],len(test_cases)))

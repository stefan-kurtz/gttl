#!/usr/bin/env python3

import argparse, sys, shutil, os
from mysubprocess import mysubprocess_expect

def files2msg_get(fail_for_protein_sequence,p):
  d = {'{}/testdata/empty.fna'.format(p)
              : ', line 1: corrupted sequence',
          '{}/testdata/fake.fna.gz'.format(p)
              : ', line 1: corrupted sequence',
          '{}/testdata/only_header.fna'.format(p)
              : ', line 2: corrupted sequence',
          '{}/testdata/only_sequence.fna'.format(p)
              : ', line 2: corrupted sequence',
          '{}/testdata/first_corrupt.fna'.format(p)
              : ', line 4: corrupted sequence',
          '{}/testdata/second_corrupt.fna'.format(p)
              : ', line 4: corrupted sequence',
          '{}/testdata/middle_of_3_err.fna'.format(p)
              : ', line 6: corrupted sequence',
          '{}/testdata/non_existing.fna'.format(p)
              : ': cannot open file',
          '{}/testdata/no_eol.fna'.format(p)
              : ', line 2: missing newline character',
         }
  if fail_for_protein_sequence:
    d['{}/testdata/protein.fsa'.format(p)] = ': can only handle DNA sequences'
  return d

def parse_arguments():
  p = argparse.ArgumentParser(description='run program for error cases')
  p.add_argument('-p','--path',type=str,default='..',
                  help='specify path of test files, default is ..')
  p.add_argument('--fail_for_protein_sequence',action='store_true',
                 default=False,help='fail for protein sequence')
  p.add_argument('-v','--valgrind',action='store_true',default=False,
                  help='run program with valgrind')
  p.add_argument('program_call',type=str,help='specify program to call')
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

files2msg = files2msg_get(args.fail_for_protein_sequence,args.path)

expected_error_code = 1
for filename, msg in files2msg.items():
  full_msg = '{}: file "{}"{}'.format(args.program_call,filename,msg)
  mysubprocess_expect('{}{} {}'.format(prefix,args.program_call,filename),
                      expected_error_code,full_msg)
sys.stderr.write('Congratulations: {} checked {} test cases for {}\n'
                  .format(sys.argv[0],len(files2msg),args.program_call))

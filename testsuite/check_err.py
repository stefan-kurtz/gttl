#!/usr/bin/env python3

import argparse, sys, shutil, os
from mysubprocess import mysubprocess_expect

def parse_arguments():
  p = argparse.ArgumentParser(description='run program for error cases')
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

files2msg \
  = {'../testdata/empty.fna' : ', line 1: corrupted sequence',
     '../testdata/only_header.fna' : ', line 2: corrupted sequence',
     '../testdata/only_sequence.fna'  : ', line 2: corrupted sequence',
     '../testdata/first_corrupt.fna'  : ', line 4: corrupted sequence',
     '../testdata/second_corrupt.fna' : ', line 4: corrupted sequence',
     '../testdata/protein.fsa' : ': can only handle DNA sequences',
     '../testdata/non_existing.fna' : ': cannot open file',
    }

expected_error_code = 1
for filename, msg in files2msg.items():
  full_msg = '{}: file "{}"{}'.format(args.program_call,filename,msg)
  mysubprocess_expect('{}{} {}'.format(prefix,args.program_call,filename),
                      expected_error_code,full_msg)
sys.stderr.write('Congratulations: {} checked {} test cases for {}\n'
                  .format(sys.argv[0],len(files2msg),args.program_call))

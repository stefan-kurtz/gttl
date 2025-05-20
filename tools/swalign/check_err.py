#!/usr/bin/env python3

import argparse, sys, shutil, os
sys.path.append('../../scripts')
from mysubprocess import mysubprocess_expect

def files2msg_get(p):
  d = {'{}/testdata/empty.fna'.format(p)
              : ', line 1: corrupted sequence',
       '{}/testdata/fake.fna.gz'.format(p)
              : ', line 1: corrupted sequence',
       '{}/testdata/only_header.fna'.format(p)
              : ', line 2: corrupted sequence',
       '{}/testdata/only_sequence.fna'.format(p)
              : ', line 1: corrupted sequence',
       '{}/testdata/first_corrupt.fna'.format(p)
              : ', line 2: corrupted sequence',
       '{}/testdata/second_corrupt.fna'.format(p)
              : ', line 4: corrupted sequence',
       '{}/testdata/middle_of_3_err.fna'.format(p)
              : ', line 4: corrupted sequence',
       '{}/testdata/non_existing.fna'.format(p)
              : ': cannot open file',
       }
  return d

def parse_arguments():
  p = argparse.ArgumentParser(description='run program for error cases')
  p.add_argument('-p','--path',type=str,default='..',
                  help='specify path of test files, default is ..')
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

files2msg = files2msg_get(args.path)

expected_error_code = 1
program_name = args.program_call.split()[0]
for filename, msg in files2msg.items():
  full_msg = '{}: file "{}"{}'.format(program_name,filename,msg)
  mysubprocess_expect('{}{} {}'.format(prefix,args.program_call,filename),
                      expected_error_code,full_msg)

err_calls_with_msgs \
  = [('{} -d testdata/gsa-seqpair.fasta -q testdata/simple_fwd.fna'.format(program_name),
      ('{}: file "testdata/gsa-seqpair.fasta": incompatible files: '
       'file "testdata/gsa-seqpair.fasta" contains protein sequences, but '
       'file "testdata/simple_fwd.fna" does not').format(program_name)),
     ('{} -d testdata/simple_fwd.fna -q testdata/gsa-seqpair.fasta'.format(program_name),
      ('{}: file "testdata/simple_fwd.fna": incompatible files: file '
       '"testdata/simple_fwd.fna" does not contain protein sequences, '
       'but file "testdata/gsa-seqpair.fasta" does').format(program_name)),
     ('{} -d testdata/seq18.fasta -g 8 1 -c 50'.format(program_name),
      ('{}: file "testdata/seq18.fasta": no Gumbel parameters for computing bits scores '
       'available for blosum62 matrix and gap '
       'parameters 8/1').format(program_name)),
     ('{} -d testdata/seq18.fasta -s unit_score_nuc'.format(program_name),
      ('{}: file "testdata/seq18.fasta": score matrix unit_score_nuc is not possible '
       'for protein sequences; the following choices are available: blosum62, '
       'unit_score_aa, unit_score_nuc, unit_score_nuc_2_2, '
       'unit_score_nuc_lower, unit_score_nuc_upper').format(program_name)),
     ('{} -d testdata/simple_fwd.fna -s unit_score_aa'.format(program_name),
      ('{}: file "testdata/simple_fwd.fna": score matrix unit_score_aa is not possible '
       'for DNA sequences; the following choices are available: blosum62, '
       'unit_score_aa, unit_score_nuc, unit_score_nuc_2_2, '
       'unit_score_nuc_lower, unit_score_nuc_upper').format(program_name)),
     ('{} -q abc'.format(program_name),
      ('{}: option -d is mandatory'.format(program_name))),
     ('{} -q'.format(program_name),
      ('{}: missing argument to option -q'.format(program_name))),
     ('{} -d'.format(program_name),
      ('{}: missing argument to option -d'.format(program_name))),
     ('{} -d testdata/seq18.fasta -g 11'.format(program_name),
      ('{}: missing arguments to option -g'.format(program_name))),
     ('{} -g 11 -d testdata/seq18.fasta'.format(program_name),
      ('{}: illegal second argument to option -g'.format(program_name))),
     ('{} -g x 11 -d testdata/seq18.fasta'.format(program_name),
      ('{}: illegal first argument to option -g'.format(program_name))),
     ('{} -t 0 -d testdata/seq18.fasta'.format(program_name),
      ('{}: argument to option -t (i.e. the number of threads) '
       'must be positive').format(program_name))] + \
     [('{} -d testdata/seq18.fasta -{}'.format(program_name,arg),
       ('{}: missing argument to option -{}'.format(program_name,arg)))
       for arg in 's'] + \
     [('{} -v {} -d testdata/seq18.fasta -s'.format(program_name,arg),
       ('{}: missing or illegal argument to option -v; '
        'argument must be either 1 or 2').format(program_name))
       for arg in ['','x','3']] + \
     [('{} -{} 0 -d testdata/seq18.fasta'.format(program_name,arg),
      ('{}: argument to option -{} must be positive'.format(program_name,arg)))
      for arg in 'bc'] + \
     [('{} -a 0 -d testdata/seq18.fasta'.format(program_name),
      (('{}: illegal argument "0" to option -a: + separated positive numbers '
        'expected').format(program_name)))] + \
     [('{} -a 2+1 -d testdata/seq18.fasta'.format(program_name),
      ('{}: illegal argument 2+1 to option -a: numbers must be strictly ordered'
       .format(program_name)))] + \
     [('{} -a 3+3 -d testdata/seq18.fasta'.format(program_name),
      ('{}: illegal argument 3+3 to option -a: numbers must be strictly ordered'
       .format(program_name)))] + \
     [('{} -a -d testdata/seq18.fasta'.format(program_name),
       '{}: missing argument to option -a'.format(program_name))] + \
     [('{} -{} -d testdata/seq18.fasta'.format(program_name,arg),
       ('{}: missing or illegal argument to option -{}'
        .format(program_name,arg)))
       for arg in 'tbc']

for call, full_msg in err_calls_with_msgs:
  mysubprocess_expect(call,expected_error_code,full_msg)

sys.stderr.write('Congratulations: {} checked {} test cases for {}\n'
                  .format(sys.argv[0],len(files2msg) + len(err_calls_with_msgs),
                          args.program_call))

#!/usr/bin/env python3

import sys, os, re, argparse
from random_sequence import RandomSequence
from popen_generator import popen_generator

def seq2fastq(sequence):
  return '\n'.join(['@',sequence,'+','~' * len(sequence)])

def nucleotide_prob(wildcard_prob):
  return (1.0 - wildcard_prob)/4

def create_sequences(tmpdirname,randseq,number,seqlen):
  if not os.path.isdir(tmpdirname):
    os.mkdir(tmpdirname)
  streams = dict()
  for format_key in 'aq':
    filename = '{}/tmp.fast{}'.format(tmpdirname,format_key)
    try:
      streams[format_key] = open(filename,'w')
    except IOError as err:
      sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
      exit(1)
  for i in range(number):
    sequence = randseq.rseq(seqlen)
    streams['a'].write('>\n{}\n'.format(sequence))
    streams['q'].write('{}\n'.format(seq2fastq(sequence)))
  for format_key in 'aq':
    streams[format_key].close()
  for format_key in 'aq':
    filename = '{}/tmp.fast{}'.format(tmpdirname,format_key)
    os.system('gzip -9 -f -k {}'.format(filename))

def run_tests(tmpdirname,program_call,max_num_threads):
  protocol = dict()
  for zipped in [True,False]:
    for format_key in 'aq':
      filename = '{}/tmp.fast{}{}'.format(tmpdirname,format_key,
                                          '.gz' if zipped else '')
      for num_threads in range(1,max_num_threads+1):
        key = '{}.{}.{}'.format(format_key,'z' if zipped else 'u',num_threads)
        protocol[key] = list()
        for line in popen_generator('{} -t {} {}'
                                    .format(program_call,num_threads,filename)):
          if not line.startswith('#'):
            protocol[key].append(line.rstrip())
        protocol[key].sort()
  keys = list(protocol.keys())
  for i in range(len(keys)-1):
    key_i = keys[i]
    key_j = keys[i+1]
    compared = '{}/{}'.format(key_i,key_j)
    if protocol[key_i] != protocol[key_j]:
      sys.stderr.write('{}: {}\n'.format(sys.argv[0],compared))
      if len(protocol[key_i]) != len(protocol[key_j]):
        sys.stderr.write('{}: protocols of different lengths {} and {}\n'
                         .format(sys.argv[0],
                                 len(protocol[key_i]),len(protocol[key_j])))
        exit(1)
      for a, b in zip(protocol[key_i],protocol[key_j]):
        if a != b:
          sys.stderr.write('{}: {} != {}\n'.format(sys.argv[0],a,b))
      exit(1)
    else:
      print('{}: identical results for "{}" and {}'.format(tmpdirname,
                                                         program_call,compared))

def parse_command_line(argv):
  default_wildcard_prob = 0.05
  default_runs = 10
  default_number_of_sequences = 1000
  default_seqlen = 300
  p = argparse.ArgumentParser(description=('run testcases of ntcard '
                                           'implementation'))
  p.add_argument('-w','--wildcard_prob',type=float,
                 default=default_wildcard_prob,
                 help=('specify probality of wildcard in random '
                       'sequences, default: {:.2f}')
                       .format(default_wildcard_prob))
  p.add_argument('-r','--runs',type=int,default=default_runs,metavar='<int>',
                 help=('specify number of testcases, default: {}'
                       .format(default_runs)))
  p.add_argument('-n','--number',type=int,default=default_number_of_sequences,
                 metavar='<int>',
                 help=('specify number of sequences per test case, default: {}'
                       .format(default_number_of_sequences)))
  p.add_argument('-l','--length',type=int,default=default_seqlen,
                 metavar='<int>',
                 help=('specify length of sequences, default: {}'
                       .format(default_seqlen)))
  p.add_argument('-t','--max_num_threads',type=int,default=3,metavar='<int>',
                 help=('specify maximum number of threads used for the test '
                       'cases, default: 3'))
  p.add_argument('program_call',type=str,default=None,
                  help='specify program call (without file argument)')
  return p.parse_args(argv)

args = parse_command_line(sys.argv[1:])
n_prob = nucleotide_prob(args.wildcard_prob)
randseq = RandomSequence('ACGTN',[n_prob] * 4 + [args.wildcard_prob])
for run in range(args.runs):
  tmpdirname = 'TMP{:02d}'.format(run)
  create_sequences(tmpdirname,randseq,args.number,args.length)
  run_tests(tmpdirname,args.program_call,args.max_num_threads)

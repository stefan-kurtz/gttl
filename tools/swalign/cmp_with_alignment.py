#!/usr/bin/env python3

import sys, re, argparse

def read_input(filename):
  try:
    stream = open(filename)
  except IOError as err:
    sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
    exit(1)
  al_lines = list()
  positions = None
  alignments = dict()
  for line in stream:
    line = line.rstrip()
    if not line.startswith('#'):
      if re.search(r'^\d',line):
        if positions is not None:
          assert len(al_lines) > 0
          alignments[positions] = al_lines.copy()
          al_lines.clear()
          positions = None
        positions = line
      else:
        al_lines.append(line)
  stream.close()
  if len(al_lines) > 0:
    alignments[positions] = al_lines.copy()
    al_lines.clear()
    positions = None
  return alignments

def parse_command_line(argv):
  p = argparse.ArgumentParser(description='compare two files with alignments')
  p.add_argument('-s','--silent',action='store_true',default=False,
                  help='do not report successfull comparison')
  p.add_argument('inputfile1',type=str,
                  help='specify first input file with alignments')
  p.add_argument('inputfile2',type=str,
                  help='specify second input file with alignments')
  return p.parse_args(argv)

args = parse_command_line(sys.argv[1:])
alignments1 = read_input(args.inputfile1)
alignments2 = read_input(args.inputfile2)
if alignments1 == alignments2:
  if not args.silent:
    print('{}: {} identical alignments'
          .format(sys.argv[0],len(alignments1)))
else:
  sys.stderr.write('{}: different alignments in {} and {}\n'
                   .format(sys.argv[0],args.inputfile1,args.inputfile2))
  if alignments1.keys() != alignments2.keys():
    sys.stderr.write('{}: different key sets in alignment in {} and {}\n'
                     .format(sys.argv[0],args.inputfile1,args.inputfile2))
    exit(1)
  for k in alignments1.keys():
    if alignments1[k] != alignments2[k]:
      sys.stderr.write('{}: different alignments for {} in {} and {}\n'
                     .format(sys.argv[0],k,args.inputfile1,args.inputfile2))
      idx = 0
      for line1,line2 in zip(alignments1[k],alignments2[k]):
        if line1 != line2:
          sys.stderr.write('first difference in line {}\n'.format(idx+1))
          sys.stderr.write('{}\n'.format(line1))
          sys.stderr.write('{}\n'.format(line2))
          exit(1)
        idx += 1
  exit(1)

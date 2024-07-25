#!/usr/bin/env python3

import sys, re, argparse
from stream_object import StreamObject

def print_sequence(seq,linelength = 70):
  for startpos in range(0,len(seq),linelength):
    print(seq[startpos:startpos+linelength])

def seq_list2sequence(seq_list):
  sequence = ' '.join(seq_list)
  return re.sub(r'\s','',sequence)

def fasta_next(filename):
  stream = StreamObject(filename)
  seq_list = list()
  header = None
  for line in stream:
    if not re.search(r'^(>|\s*$|\s*#)',line):
      seq_list.append(line)
    elif len(line) > 0 and line[0] == '>':
      if header is not None:
        yield header, seq_list2sequence(seq_list)
        seq_list.clear()
      header = line[1:]
  if header is not None:
    yield header, seq_list2sequence(seq_list)

def fastq_next(filename):
  stream = StreamObject(filename)
  state = 0
  header = None
  for line in stream:
    line = line.rstrip()
    if state == 0:
      assert len(line) > 0 and line[0] == '@'
      header = line[1:]
      state += 1
    elif state == 1:
      yield header, line
      state += 1
    elif state == 2:
      state += 1
    else:
      state = 0

def reader_get(filename):
  if re.search(r'\.fq$',filename) or re.search(r'\.fastq$',filename):
    return fastq_next
  return fasta_next

def parse_command_line(argv):
  p = argparse.ArgumentParser(description='parse fasta or fastq file')
  p.add_argument('-w','--width',type=int,default=64,metavar='<int>',
                 help='specify width of lines in sequence output')
  p.add_argument('inputfiles',nargs='+',
                 help='specify input files, - means stdin')
  return p.parse_args(argv)

if __name__ == '__main__':
  args = parse_command_line(sys.argv[1:])
  for filename in args.inputfiles:
    reader = reader_get(filename)
    for header, sequence in reader(filename):
      print('>{}'.format(header))
      print_sequence(sequence, args.width)

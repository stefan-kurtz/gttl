#!/usr/bin/env python3

import sys, re, argparse

def print_sequence(seq,linelength = 70):
  for startpos in range(0,len(seq),linelength):
    print(seq[startpos:startpos+linelength])

def seq_list2sequence(seq_list):
  sequence = ' '.join(seq_list)
  return re.sub('\s','',sequence)

def stream_get(filename):
  if filename == '-':
    return sys.stdin
  try:
    stream = open(filename,'r')
  except IOError as err:
    sys.stderr.write('{}: {}\n'
                      .format(sys.argv[0],err))
    exit(1)
  return stream

def fasta_next(filename):
  stream = stream_get(filename)
  seq_list = list()
  header = None
  for line in stream:
    line = line.rstrip()
    if not re.search(r'^(>|\s*$|\s*#)',line):
      seq_list.append(line.rstrip())
    elif len(line) > 0 and line[0] == '>':
      if header is not None:
        yield header, seq_list2sequence(seq_list)
        seq_list.clear()
      header = line[1:]
  if header is not None:
    yield header, seq_list2sequence(seq_list)
  if filename != '-':
    stream.close()

def fastq_next(filename):
  stream = stream_get(filename)
  seq_list = list()
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
  if filename != '-':
    stream.close()

def parse_arguments():
  p = argparse.ArgumentParser(description='parse fasta or fastq file')
  p.add_argument('-w','--width',type=int,default=64,metavar='<int>',
                 help='specify width of lines in sequence output')
  p.add_argument('inputfiles',nargs='+',
                 help='specify input files, - means stdin')
  return p.parse_args()

if __name__ == '__main__':
  args = parse_arguments()
  for filename in args.inputfiles:
    if re.search(r'\.fq$',filename) or re.search(r'\.fastq$',filename):
      reader = fastq_next
    else:
      reader = fasta_next
    for header, sequence in reader(filename):
      print('>{}'.format(header))
      print_sequence(sequence, args.width)

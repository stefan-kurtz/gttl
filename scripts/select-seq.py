#!/usr/bin/env python3

import sys, re, argparse, random
import fastaIterator
from print_sequence import print_sequence

def keyset_from_file(key):
  assert key
  try:
    stream = open(key)
  except IOError as err:
    return set([key])
  key_set = set()
  for line in stream:
    key_set.add(line.strip())
  stream.close()
  return key_set

def seqnum_extract_map_get(args_number):
  sequence_extract_map = dict()
  for line in keyset_from_file(args_number):
    assert len(line) > 0
    arr = line.split()
    assert len(arr) > 0
    sequence_number = int(arr[0])
    if sequence_number in sequence_extract_map:
      sys.stderr.write(('{}: entry for sequence number {} in file {} already '
                        'exists\n').format(sys.argv[0],sequence_number,
                                           args_number))
      exit(1)
    this_substring_length_range = None
    if len(arr) == 3:
      this_substring_length_range = (int(arr[1]),int(arr[2]))
    elif len(arr) == 2:
      start_pos = int(arr[1])
      this_substring_length_range = (start_pos,None)
    elif len(arr) == 1:
      this_substring_length_range = (0,None)
    else:
      sys.stderr.write(('{}: option -n must be a string of 1-3 int arguments, '
                        'or if the argument is a filename all line must '
                        'consist of 1-3 int arguments\n')
                       .format(sys.argv[0]))
      exit(1)
    sequence_extract_map[sequence_number] = this_substring_length_range
  return sequence_extract_map

def check_positive(value):
  try:
    ivalue = int(value)
  except ValueError as err:
    raise argparse.ArgumentTypeError('{} is an int value'.format(value))
  if ivalue <= 0:
    raise argparse.ArgumentTypeError('{} is an invalid positive int value'
                                      .format(value))
  return ivalue

def number_of_sequences(inputfiles):
  count = 0
  for inputfile in inputfiles:
    reader = fastaIterator.reader_get(inputfile)
    for header,sequence in reader(inputfile):
      count += 1
  return count

def parse_command_line():
  p = argparse.ArgumentParser()
  group = p.add_mutually_exclusive_group(required=True)
  group.add_argument('-f','--first',type=check_positive,
                     metavar='<n>',default=None,
                     help='select the first <n> sequences')
  group.add_argument('-l','--last',type=check_positive,
                     metavar='<n>',default=None,
                     help='select the last <n> sequences')
  group.add_argument('-n','--number',type=str,
                     metavar='<n> [startpos] [len]',default=None,
                     help=('select sequence with number <n> (count from 0), '
                           'optionally specify the start position and '
                           '(optionally) the length of the sequence to be '
                           'extracted; if the length is omitted, '
                           'the suffix of the sequence is taken; specify '
                           'all numbers in one string, like \'23 42 20\' '
                           'to extract from sequence number 23 the substring '
                           'of length 20 begininng at position 42; like '
                           '\'42 300\' to extract from sequence 42 the suffix '
                           'beginning at position 300; all positions are '
                           '0-based'))
  group.add_argument('-r','--random_sample',type=int,
                     metavar='<n>',default=None,
                     help=('output a random sample of the specified number of '
                           'sequences'))
  group.add_argument('-p','--pattern',type=str,metavar='<p>',default=None,
                     help=('select sequences matching pattern <p> in header '
                           '(default) or sequence (if option --s is used) '
                           'if <p> is a path to a file, then patterns are '
                           'taken from file, one per line'))
  group.add_argument('--lrange',type=int,nargs=2,action='store',
                     metavar='<N>',
                     default=None,
                     help=('select sequences of given minimum and maximum '
                           'length'))
  searchwhere = p.add_mutually_exclusive_group(required=False)
  searchwhere.add_argument('-i','--identifier',type=str,default=None,
                           choices=['yes','no'],
                           help=('search pattern in ID of header, i.e. '
                                 'leftmost string delimited by | on both '
                                 'sides, '))
  searchwhere.add_argument('-s','--seq',action='store_true',default=False,
                           help='search pattern in sequence rather than header')
  p.add_argument('--substring_length',type=str,default=None,
                 help=('only output substring of some length in the given '
                       'range selected at a random position from each selected '
                       'sequence'))
  p.add_argument('-c','--case_insensitive',action='store_true',default=False,
                 help=('match case insensitive'))
  p.add_argument('-w','--width',type=int,
                 default=60,
                 help=('specify width of line of formatted sequence, '
                       'default=60'))
  p.add_argument('--fastq_output',action='store_true',default=False,
                  help=('output sequences in fastq format with fake quality '
                        'lines'))
  p.add_argument('inputfiles',nargs='+',
                  help='specify input files, - means stdin')
  return p.parse_args()

def this_print_sequence(seq,linelength,stream = sys.stdout):
  if linelength is None:
    linelength = len(seq)
  for startpos in range(0,len(seq),linelength):
    stream.write('{}\n'.format(seq[startpos:startpos+linelength]))

def show_the_sequence(header,sequence,substring_length_range,width,
                      fastq_output,no_random=False):
  if fastq_output:
    header_mark = '@'
    width = None
  else:
    header_mark = '>'
  if substring_length_range:
    if len(sequence) >= substring_length_range[1]:
      if no_random:
        r_start = substring_length_range[0]
        r_len = substring_length_range[1]
      else:
        r_len = random.randint(substring_length_range[0],
                               substring_length_range[1])
        r_start = random.randint(0,len(sequence) - r_len)
      print('{}{} [{},{}]'.format(header_mark,header,r_start,r_start+r_len-1))
      sequence = re.sub(r'[XUBZJO]+','',sequence)
      this_print_sequence(sequence[r_start:r_start+r_len],width)
      if fastq_output:
        print('+')
        this_print_sequence('E' * r_len,width)
  else:
    print('{}{}'.format(header_mark,header))
    sequence = re.sub(r'[XUBZJO]+','',sequence)
    this_print_sequence(sequence,width)
    if fastq_output:
      print('+')
      this_print_sequence('E' * len(sequence),width)

args = parse_command_line()
last_queue = list()
if args.pattern:
  pattern_set = keyset_from_file(args.pattern)
else:
  pattern_set = None

seq_count = 0
stop = False
sample_seq_nums = None
substring_length_range = None
sequence_extract_map = None
if args.substring_length:
  mo = re.search(r'(\d+)-(\d+)',args.substring_length)
  if mo:
    substring_length_range = (int(mo.group(1)),int(mo.group(2)))
  else:
    sys.stderr.write('{}: illegal length range specified with option '
                     ' --substring_length\n'
                     .format(sys.argv[0]))
    exit(1)
if args.random_sample:
  num_seq = number_of_sequences(args.inputfiles)
  if args.random_sample > num_seq:
     sys.stderr.write(('{}: cannot sample {} sequences from a file with {} '
                       'sequences\n').format(sys.argv[0],args.random_sample,
                                             num_seq))
     exit(1)
  seq_nums = list(range(num_seq))
  sampled_seq_nums = sorted(random.sample(seq_nums,args.random_sample))
for inputfile in args.inputfiles:
  if stop:
    break
  reader = fastaIterator.reader_get(inputfile)
  for header,sequence in reader(inputfile):
    okay = False
    if not (args.random_sample is None):
      if len(sampled_seq_nums) == 0:
        stop = True
        break
      if seq_count == sampled_seq_nums[0]:
        show_the_sequence(header,sequence,substring_length_range,args.width,
                          args.fastq_output)
        sampled_seq_nums.pop(0)
      else:
        assert seq_count < sampled_seq_nums[0]
    if not (args.first is None):
      if seq_count >= args.first:
        stop = True
        break
      else:
        show_the_sequence(header,sequence,substring_length_range,args.width,
                          args.fastq_output)
    if not (args.last is None):
      if len(last_queue) == args.last:
        last_queue.pop(0)
      assert len(last_queue) < args.last
      last_queue.append([header,sequence])
    elif not (pattern_set is None):
      if args.identifier:
        idm = re.search(r'\|([^\|]+)\|',header)
        if not idm:
          sys.stderr.write('{}: cannot find identifier in {}\n'
                            .format(sys.argv[0],header))
          exit(1)
        ident = idm.group(1)
        if (args.identifier == 'yes' and ident in pattern_set) or \
           (args.identifier == 'no' and not ident in pattern_set):
          show_the_sequence(header,sequence,substring_length_range,args.width,
                            args.fastq_output)
      else:
        for pattern in pattern_set:
          if args.seq:
            where_to_look = sequence
          else:
            where_to_look = header
          if (args.case_insensitive and \
              re.search(r'{}'.format(pattern),where_to_look,flags=re.I)) or \
             (not args.case_insensitive and \
              re.search(r'{}'.format(pattern),where_to_look)):
            show_the_sequence(header,sequence,substring_length_range,args.width,
                              args.fastq_output)
    elif not (args.number is None):
      if sequence_extract_map is None:
        sequence_extract_map = seqnum_extract_map_get(args.number)
      if seq_count in sequence_extract_map:
        this_substring_length_range = sequence_extract_map[seq_count]
        if this_substring_length_range[1] is None:
          this_substring_length_range = (this_substring_length_range[0],
                                         len(sequence) - 
                                         this_substring_length_range[0])
        no_random = True
        show_the_sequence(header,sequence,this_substring_length_range,
                          args.width,args.fastq_output,no_random)
    elif not (args.lrange is None):
      minlen = args.lrange[0]
      maxlen = args.lrange[1]
      seqlen = len(sequence)
      if seqlen >= minlen and seqlen <= maxlen:
        show_the_sequence(header,sequence,substring_length_range,args.width,
                          args.fastq_output)
    seq_count += 1

if not (args.last is None):
  for header, sequence in last_queue:
    show_the_sequence(header,sequence,substring_length_range,args.width,
                      args.fastq_output)

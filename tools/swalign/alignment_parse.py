#!/usr/bin/env python3

import sys, re, argparse

class AlignmentLine:
  def __init__(self,mo):
    self.startpos = int(mo.group(1))
    self.aligned_seq = mo.group(2)
    self.endpos = int(mo.group(3))

def num_aligned(s):
 return len(s) - s.count('-')

class AlignmentUnit:
  def __init__(self,coords):
    self.coords = coords
    self.subject_lines = list()
    self.middle_lines = list()
    self.query_lines = list()
  def append_subject_line(self,subject_line):
    self.subject_lines.append(subject_line)
  def append_middle_line(self,middle_line):
    self.middle_lines.append(middle_line)
  def append_query_line(self,query_line):
    self.query_lines.append(query_line)
  def coords_get(self):
    return self.coords
  def subject_aligned(self):
    return ''.join([s.aligned_seq for s in self.subject_lines])
  def query_aligned(self):
    return ''.join([q.aligned_seq for q in self.query_lines])
  def middle_line(self):
    return ''.join(self.middle_lines)
  def verify(self):
    s_start = int(self.coords['s. start'])
    s_len = int(self.coords['s. len'])
    s_end = s_start + s_len - 1
    q_start = int(self.coords['q. start'])
    q_len = int(self.coords['q. len'])
    q_end = q_start + q_len - 1
    assert s_start == self.subject_lines[0].startpos, \
           '{} != {}'.format(s_start,self.subject_lines[0].startpos)
    assert s_end == self.subject_lines[-1].endpos, \
           '{} != {}'.format(s_end,self.subject_lines[-1].endpos)
    assert q_start == self.query_lines[0].startpos, \
           '{} != {}'.format(q_start,self.query_lines[0].startpos)
    assert q_end == self.query_lines[-1].endpos, \
           '{} != {}'.format(q_end,self.query_lines[-1].endpos)
    assert len(self.subject_aligned()) == len(self.query_aligned()), \
           '{} != {}'.format(len(self.subject_aligned()) == \
                             len(self.query_aligned()))
    assert len(self.subject_aligned()) == len(self.middle_line()), \
           '{} != {}'.format(len(self.subject_aligned()) == \
                             len(self.middle_line()))
    for items in [self.subject_lines,self.query_lines]:
      previous = None
      for p_val in items:
        assert previous is None or previous.endpos + 1 == p_val.startpos
        assert p_val.startpos + num_aligned(p_val.aligned_seq) == p_val.endpos+1
        previous = p_val

def split_num_seq(line):
  mo = re.search(r'^\s*(\d+)\s*(\S+)\s*(\d+)\s*$',line)
  if not mo:
    sys.stderr.write(('{}: cannot parse {} as number enclosed aligned '
                      'sequence\n').format(sys.argv[0],line))
    exit(1)
  return AlignmentLine(mo)

def determine_indent(key,line):
  mo = re.search(r'^({}\s+\d+\s*)\S'.format(key),line)
  if not mo:
    sys.stderr.write('{}: cannot parse number prefix from {}'
                     .format(sys.argv[0],line))
    exit(1)
  return len(mo.group(1))

def enum_alignments(inputfile):
  if inputfile != '-':
    try:
      stream = open(inputfile)
    except IOError as err:
      sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
      exit(1)
  else:
    stream = sys.stdin
  fields = None
  options_prefix = '# Options: '
  alignment = None
  fields = None
  options = None
  subject_indent = None
  for line in stream:
    line = line.rstrip('\n')
    field_pattern = '# Fields: '
    if line.startswith(field_pattern):
      if fields is None:
        fields = line[len(field_pattern):].split(', ')
    elif line.startswith(options_prefix):
      if options is None:
        options = line[len(options_prefix):]
    elif not line.startswith('#'):
      if len(line) > 0:
        if line.startswith('Sbjct'):
          subject_indent = determine_indent('Sbjct',line)
          alignment.append_subject_line(split_num_seq(line[6:]))
        elif line.startswith('Query'):
          query_indent = determine_indent('Query',line)
          assert subject_indent == query_indent
          alignment.append_query_line(split_num_seq(line[6:]))
        elif line[0] == ' ':
          alignment.append_middle_line(line[subject_indent:])
        else:
          if alignment is not None:
            yield alignment
          arr = line.split('\t')
          assert len(arr) >= 2, 'arr={}'.format(arr)
          alignment = AlignmentUnit({f:v for f,v in zip(fields,arr)})
  if alignment is not None:
    yield alignment
  if inputfile != '-':
    stream.close()

def parse_command_line(argv):
  p = argparse.ArgumentParser(description=('collect information about '
                                           'formatted alignments'))
  p.add_argument('inputfile',type=str,
                  help='specify input file, means stdin')
  return p.parse_args(argv)

args = parse_command_line(sys.argv[1:])
for idx, alignment in enumerate(enum_alignments(args.inputfile)):
  alignment.verify()

import sys, re, os
from stream_object import StreamObject
from intervals_compare import coords_contained_in, coords_overlap_size

class SEmatch:
  def __init__(self,matchinput,force_some = True):
    self._matchinput = matchinput
    self._force_some = force_some
    self._options = None
    self._fields = None
    self._runtime = None
    self._spacepeak = None
  def each(self):
    stream_object = StreamObject(self._matchinput)
    for line in stream_object:
      if line.startswith('#'):
        mo = re.search(r'^# Options:(.*)$',line)
        if mo is not None:
          self._options = mo.group(1)
        else:
          mo = re.search(r'^# Fields:\s(.*)$',line)
          if mo is not None:
            fields = re.sub(r'\.','_',mo.group(1)).split(', ')
            self._fields = list()
            for f in fields:
              f1 = re.sub(r'\s','_',f)
              f2 = re.sub(r'__','_',f1)
              f3 = re.sub(r'%_','',f2)
              self._fields.append(f3)
          else:
            mo = re.search(r'^# TIME.*\s(\S+)$',line)
            if mo:
              self._runtime = float(mo.group(1))
            else:
              mo = re.search(r'^# combined space peak in megabytes:\s(\S+)$',
                             line)
              if mo:
                self._spacepeak = float(mo.group(1))
      else:
        if self._fields is None:
          sys.stderr.write('{}: {}: fields: not defined\n'
                            .format(sys.argv[0],self._matchinput))
          exit(1)
        val_dict = dict()
        line = re.sub(r'^\s+','',line)
        for idx, value in enumerate(line.split()):
          if re.search(r'^\d+\.\d+$',value):
            val_dict[self._fields[idx]] = float(value)
          elif re.search(r'^\d+$',value):
            val_dict[self._fields[idx]] = int(value)
          else:
            val_dict[self._fields[idx]] = value
        if 's_start' in val_dict:
          if 's_len' in val_dict:
            assert isinstance(val_dict['s_start'],int)
            assert isinstance(val_dict['s_len'],int)
            val_dict['s_end'] = val_dict['s_start'] + val_dict['s_len'] - 1
          elif 's_end' in val_dict:
            if val_dict['s_start'] < val_dict['s_end']:
              val_dict['s_len'] = val_dict['s_end'] - val_dict['s_start'] + 1
            else:
              val_dict['s_len'] = val_dict['s_start'] - val_dict['s_end'] + 1
          else:
            sys.stderr.write(('{}: either length of match on subject or end '
                              'position must be given\n').format(sys.argv[0]))
            exit(1)
        else:
          if self._force_some:
            sys.stderr.write('{}: start of match on subject must be given'
                              .format(sys.argv[0]))
            exit(1)
        if 'q_start' in val_dict:
          if 'q_len' in val_dict:
            val_dict['q_end'] = val_dict['q_start'] + val_dict['q_len'] - 1
          elif 'q_end' in val_dict:
            if val_dict['q_start'] < val_dict['q_end']:
              val_dict['q_len'] = val_dict['q_end'] - val_dict['q_start'] + 1
            else:
              val_dict['q_len'] = val_dict['q_start'] - val_dict['q_end'] + 1
          else:
            val_dict['q_len'] = val_dict['s_len']
        else:
          if self._force_some:
            sys.stderr.write('{}: start of match on query must be given\n'
                              .format(sys.argv[0]))
            exit(1)
        val_dict['origline'] = line.rstrip()
        yield val_dict
  def spacepeak_get(self):
    return self._spacepeak
  def runtime_get(self):
    return self._runtime
  def options_get(self):
    return self._options
  def fields_get(self):
    return self._fields

def match_is_identical(m0,m1):
  if (m0['s_seqnum'] != m1['s_seqnum']) or (m0['q_seqnum'] != m1['q_seqnum']):
    sys.stderr.write('{}: expect same sequence numbers\n'.format(sys.argv[0]))
    exit(1)
  if [m0['s_start'],
      m0['s_end'],
      m0['q_start'],
      m0['q_end']] == [m1['s_start'],m1['s_end'],m1['q_start'],m1['q_end']]:
    return True
  return False

def coords_oversize(start0,end0,start1,end1):
  assert coords_contained_in(start0,end0,start1,end1)
  return (start0 - start1) + (end1 - end0)

def match_proper_contained_in(m0,m1):
  if m0['s_seqnum'] != m1['s_seqnum'] or m0['q_seqnum'] != m1['q_seqnum']:
    sys.stderr.write('{}: expect same sequence numbers\n'.format(sys.argv[0]))
    exit(1)
  if coords_contained_in(m0['s_start'],m0['s_end'],m1['s_start'],m1['s_end']) \
     and \
     coords_contained_in(m0['q_start'],m0['q_end'],m1['q_start'],m1['q_end']) \
     and \
     not match_is_identical(m0,m1):
    return True
  return False

def match_oversize(m0,m1):
  if m0['s_seqnum'] != m1['s_seqnum'] or m0['q_seqnum'] != m1['q_seqnum']:
    sys.stderr.write('{}: expect same sequence numbers\n'.format(sys.argv[0]))
    exit(1)
  return (coords_oversize(m0['s_start'],m0['s_end'],m1['s_start'],m1['s_end']) +
          coords_oversize(m0['q_start'],m0['q_end'],m1['q_start'],m1['q_end']))

def matches_overlap(m0,m1):
  ovl = 0
  if m0['s_start'] <= m1['s_start']:
    ovl +=  coords_overlap_size(m0['s_start'],m0['s_end'],
                                m1['s_start'],m1['s_end'])
  else:
    ovl +=  coords_overlap_size(m1['s_start'],m1['s_end'],
                                m0['s_start'],m0['s_end'])
  if m0['q_start'] <= m1['q_start']:
    ovl +=  coords_overlap_size(m0['q_start'],m0['q_end'],
                                m1['q_start'],m1['q_end'])
  else:
    ovl +=  coords_overlap_size(m1['q_start'],m1['q_end'],
                                m0['q_start'],m0['q_end'])
  len_val = (m0['s_len'] + m0['q_len'])/2 + (m1['s_len'] + m1['q_len'])/2
  return round(float(ovl)/float(len_val))

def match_gap(end0,start1):
  assert end0 < start1, ('{}: end0={} >= {}=start1\n'
                         .format(sys.argv[0],end0,start1))
  return start1 - (end0 + 1)

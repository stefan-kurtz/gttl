#!/usr/bin/env python3

import sys, re

def file2keyhash(filename):
  if filename != '-':
    try:
      stream = open(filename)
    except IOError as err:
      sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
      exit(1)
  else:
    stream = sys.stdin
  keyhash = dict()
  seqpairs = list()
  for line in stream:
    if not re.search(r'^#',line):
      ll = line.split()
      if len(ll) != 8:
        sys.stderr.write('{}: file {}: expect 8 columns per line\n'
                          .format(sys.argv[0],filename))
        exit(1)
      sseq = int(ll[0])
      qseq = int(ll[1])
      key = int(ll[3]) + int(ll[5])
      if not (sseq in  keyhash):
        keyhash[sseq] = dict()
      if qseq in keyhash[sseq]:
        sys.stderr.write('{}: sequence pair {} {} occurs twice'
                          .format(sys.argv[0],sseq,qseq))
        exit(1)
      seqpairs.append((sseq,qseq))
      keyhash[sseq][qseq] = key
  if filename != '-':
    stream.close()
  return keyhash, sorted(seqpairs)

if len(sys.argv) != 3:
  sys.stderr.write('Usage: {} <matchfile1> <matchfile2>\n'.format(sys.argv[0]))
  exit(1)

keyhash1, seqpairs1 = file2keyhash(sys.argv[1])
keyhash2, seqpairs2 = file2keyhash(sys.argv[2])

if seqpairs1 != seqpairs2:
  sys.stderr.write('{}: different seqpair lists\n'.format(sys.argv[0]))
  exit(1)

firstlarger = 0
secondlarger = 0
for seqpair in  seqpairs1:
  sseq, qseq = seqpair
  v1 = keyhash1[sseq][qseq]
  v2 = keyhash2[sseq][qseq]
  if v1 > v2:
    sys.stderr.write('{}: for key {}/{}: v1 = {} > {} = v2\n'
                      .format(sys.argv[0],sseq,qseq,v1,v2))
    firstlarger += 1
  elif v1 < v2:
    sys.stderr.write('{}: for key {}/{}: v1 = {} < {} = v2\n'
                      .format(sys.argv[0],sseq,qseq,v1,v2))
    secondlarger += 1

if firstlarger > 0 or secondlarger > 0:
  if firstlarger > 0:
    sys.stderr.write("firstlarger={}\n".format(firstlarger))
  if secondlarger > 0:
    sys.stderr.write("secondlarger={}\n".format(secondlarger))
  exit(1)

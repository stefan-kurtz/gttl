#!/usr/bin/env python3

import sys, re, argparse, os

def filepath2tag(filepath):
  filename = os.path.basename(filepath)
  return re.sub(r'\.','_',filename.upper())

def parse_command_line(argv):
  p = argparse.ArgumentParser(description=('check for all given input files, '
                                           'that ifndef/define in .hpp files '
                                           'is consistent with filename'))
  p.add_argument('inputfiles',nargs='+', help='specify input files')
  return p.parse_args(argv)

args = parse_command_line(sys.argv[1:])

for filepath in args.inputfiles:
  try:
    stream = open(filepath)
  except IOError as err:
    sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
    exit(1)
  expected_tag = filepath2tag(filepath)
  found_ifndef = False
  found_define = False
  for line in stream:
    line = line.rstrip()
    mo = re.search(r'^#ifndef\s+(.*)$',line)
    if mo:
      found_ifndef = True
      tag = mo.group(1)
      if tag != expected_tag:
        sys.stderr.write('{}: file {}, illegal line\n{}\n'
                         .format(sys.argv[0],filepath,line))
        sys.stderr.write('must be\n#ifndef {}\n'.format(expected_tag))
        exit(1)
      found_ifndef = True
    else:
      mo = re.search(r'^#define\s+(.*)$',line)
      if mo:
        if not found_ifndef:
          sys.stderr.write('{}: file {}: illegal line {}\n'
                           .format(sys.argv[0],filepath,line))
          sys.stderr.write('missing line\n#ifndef {}\n'.format(expected_tag))
          exit(1)
        tag = mo.group(1)
        if tag != expected_tag:
          sys.stderr.write('{}: file {}: illegal line {}\n'
                           .format(sys.argv[0],filepath,line))
          sys.stderr.write('must be\n#ifndef {}\n'.format(expected_tag))
          exit(1)
        found_define = True
        break
  if not found_ifndef:
    sys.stderr.write('{}: file {}: missing line #ifndef {}\n'
                     .format(sys.argv[0],filepath,expected_tag))
    exit(1)
  if not found_define:
    sys.stderr.write('{}: file {}: missing line #define {}\n'
                     .format(sys.argv[0],filepath,expected_tag))
    exit(1)

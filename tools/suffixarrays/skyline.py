#!/usr/bin/env python3

import sys, argparse

def skyline(i,reverse = False):
  prev_s = str()
  char_map = [chr(ord('a') + j) for j in range(0,i)]
  if reverse:
    char_map = list(reversed(char_map))
  for j in range(1,i+1):
    cc = char_map[j-1]
    s = "{}{}{}".format(prev_s,cc,prev_s)
    prev_s = s
  return prev_s

def find_level(alphasize,j):
  assert alphasize > 0
  mid = 2**(alphasize-1) - 1
  level = 0
  while mid > 0 and j != mid:
    if j > mid:
      j -= (mid+1)
    level += 1
    mid //= 2
  return level

def skyline_at(alphasize,s):
  assert len(s) == 2**alphasize - 1
  char_map = [chr(ord('a') + j) for j in range(0,alphasize)]
  for j, cc in enumerate(s):
    l = find_level(alphasize,j)
    assert l < alphasize and char_map[alphasize-1-l] == cc

def parse_arguments():
  p = argparse.ArgumentParser(description='output skyline string')
  p.add_argument('-f','--fasta',action='store_true',default=False,
                 help='show output in fasta format')
  p.add_argument('-n','--no_nl',action='store_true',default=False,
                 help='do not append newline at end of sequence')
  p.add_argument('-s','--sentinel',action='store_true',default=False,
                 help=('add sentinel larger than any other character '
                       ' at end of sequence'))
  p.add_argument('--pentinel',action='store_true',default=False,
                 help=('add sentinel smaller than any other character '
                       ' at end of sequence'))
  p.add_argument('-r','--reverse',action='store_true',default=False,
                 help='reverse character order')
  p.add_argument('--skyline_at',action='store_true',default=False,
                 help=('call skyline_at function to verify the mapping from '
                       'indexes of the skyline string to their contents'))
  p.add_argument('alphasize',type=int,
                  help='specify alphabet size of skyline string')
  return p.parse_args()

args = parse_arguments()

sstring = skyline(args.alphasize,args.reverse)
if args.skyline_at:
  skyline_at(args.alphasize,sstring)
if args.fasta:
  print('>\n{}'.format(sstring),end='')
else:
  print(sstring,end='')

if args.sentinel:
  print(chr(127),end='')

if args.pentinel:
  print('Z',end='')

if not args.no_nl:
  print("")

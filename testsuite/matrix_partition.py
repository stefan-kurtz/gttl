#!/usr/bin/env python3

import sys, argparse, functools
import numpy as np

def split_itv(a,b):
  h = b//2 + (b % 2)
  return [(a,h),(a+h,b-h)]

def partition_two(cutlen,m,n):
  stack = [((0,m),(0,n))]
  while len(stack) > 0:
    ((i,j),(k,l)) = stack.pop()
    if j <= cutlen and l <= cutlen:
      yield ((i,j),(k,l))
    else:
      if j < l:
        for itv in split_itv(k,l):
          stack.append(((i,j),itv))
      else:
        for itv in split_itv(i,j):
          stack.append((itv,(k,l)))

def partition_interval(cutlen,m):
  for i in range(0,m,cutlen):
    if i + cutlen <= m:
      yield (i,cutlen)
    else:
      if m > i:
        yield (i,m - i)

def partition_self(cutlen,m):
  ps = list(partition_interval(cutlen,m))
  for p in ps:
    yield (p,None)
  for i in range(len(ps)):
    for j in range(i+1,len(ps)):
      yield (ps[i],ps[j])

def antidiagonal(itv):
  if itv[1] is None:
    return 2 * itv[0][0]
  else:
    return itv[0][0] + itv[1][0]

def antidiagonal_cmp(l,r):
  if antidiagonal(l) < antidiagonal(r):
    return -1
  if antidiagonal(l) > antidiagonal(r):
    return +1
  if l[0][0] < r[0][0]:
    return -1
  if l[0][0] > r[0][0]:
    return 1
  assert false

def partition(cutlen,m,n):
  l = list()
  if n == 0:
    for v in partition_self(cutlen,m):
      l.append(v)
  else:
    for v in partition_two(cutlen,m,n):
      l.append(v)
  return sorted(l,key=functools.cmp_to_key(antidiagonal_cmp))

def pairs_count(covered,p):
  s0 = p[0][0]
  l0 = p[0][1]
  end0 = s0 + l0
  if p[1] is None:
    if covered is not None:
      for i in range(s0,end0-1):
        for j in range(i+1,end0):
          assert covered[i][j] == 0
          covered[i][j] = 1
    return (l0 * (l0-1))//2
  else:
    s1 = p[1][0]
    l1 = p[1][1]
    end1 = s1 + l1
    if covered is not None:
      for i in range(s0,end0):
        for j in range(s1,end1):
          assert covered[i][j] == 0
          covered[i][j] = 1
    return l0 * l1

def parse_arguments(argv):
  p = argparse.ArgumentParser(description='divide matrix into submatrices')
  p.add_argument('-d','--debug',action='store_true',default=False,
                  help=('keep track of covered intervals (very slow) for '
                        'debugging'))
  p.add_argument('cutlen',type=int,default=None,help='specify cutlen')
  p.add_argument('rows',type=int,default=None,help='specify number of rows')
  p.add_argument('cols',type=int,default=None,help='specify number of columns')
  return p.parse_args(argv)

args = parse_arguments(sys.argv[1:])

collect = set()
num_pairs = 0
covered = None
if args.debug:
  if args.cols == 0:
    covered = np.zeros((args.rows,args.rows),dtype=np.uint8)
  else:
    covered = np.zeros((args.rows,args.cols),dtype=np.uint8)

if args.cols == 0:
  expected = (args.rows * (args.rows-1))//2
else:
  expected = args.rows * args.cols

for itv in partition(args.cutlen,args.rows,args.cols):
  assert not itv in collect
  collect.add(itv)
  pc = pairs_count(covered,itv)
  num_pairs += pc
  print('{}\t{}\t{}'.format(itv,antidiagonal(itv),pc))

print('# number of parts\t{}'.format(len(collect)))
print('# number of pairs\t{}'.format(num_pairs))
assert not args.debug or num_pairs == np.sum(covered)
assert num_pairs == expected, ('num_pairs = {} != {} = expected'
                               .format(num_pairs,expected))

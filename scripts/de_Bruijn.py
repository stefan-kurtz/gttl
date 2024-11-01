#!/usr/bin/env python3

import argparse

'''
A de Bruijn sequence of order q over an alphabet A is a cyclic string in
which every possible q-gram over A occurs exactly once as a substring.
The following Python code calculates a de Bruijn sequence for a given
alphabet and q, based on an algorithm from Frank Ruskey's Combinatorial
Generation (see
https://git.sagemath.org/sage.git/tree/src/sage/combinat/debruijn_sequence.pyx)
'''

def de_bruijn(alphabet: str, q: int) -> str:
    alpha_size = len(alphabet)
    a = [0] * alpha_size * q
    sequence = []

    def dbseq(t, p):
        if t > q:
            if q % p == 0:
                sequence.extend(a[1 : p + 1])
        else:
            a[t] = a[t - p]
            dbseq(t + 1, p)
            for j in range(a[t - p] + 1, alpha_size):
                a[t] = j
                dbseq(t + 1, t)

    dbseq(1, 1)
    return ''.join(alphabet[i] for i in sequence)

def parse_arguments():
  p = argparse.ArgumentParser(description=('generate de Bruijn sequence with '
                                           'respect to q and given alphabet'))
  p.add_argument('-q','--qgram_length',type=int,default=3,metavar='<int>',
                  help='specify parameter q')
  alphabet_group = p.add_mutually_exclusive_group(required=True)
  p.add_argument('-t','--transform',action='store_true',default=False,
                  help=('use integer code from 0 to alpha_size-1 for each '
                        'character'))
  alphabet_group.add_argument('-a','--alphabet',type=str,metavar='<characters>',
                              help='specify alphabet')
  alphabet_group.add_argument('--dna',action='store_true',default=False,
                              help='use DNA alphabet')
  alphabet_group.add_argument('--protein',action='store_true',default=False,
                              help='use Protein alphabet')
  p.add_argument('--fasta',action='store_true',default=False,
                 help='output sequences in fasta format')
  return p.parse_args()

args = parse_arguments()
if args.dna:
  args.alphabet = 'ACGT'
elif args.protein:
  args.alphabet = 'ADEFGHIKLMNPQRSTVWY'
if args.transform:
  args.alphabet = ''.join(map(chr,range(len(args.alphabet))))
seq = de_bruijn(args.alphabet,args.qgram_length)
if args.fasta:
  print('>de Bruijn sequence for alphabet={}, q={}'.format(args.alphabet,args.qgram_length))
  width = 60
  for start in range(0,len(seq),width):
    print('{}'.format(seq[start:start + width]))
else:
  print(seq,end='')

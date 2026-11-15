#!/usr/bin/env python3

'''
A de Bruijn sequence of order q over an alphabet A is a cyclic string in
which every possible q-gram over A occurs exactly once as a substring.
The following Python code calculates a de Bruijn sequence for a given
alphabet and q, based on an algorithm from Frank Ruskey's Combinatorial
Generation (see
https://git.sagemath.org/sage.git/tree/src/sage/combinat/debruijn_sequence.pyx)
The cyclic sequence is extended by its first q-1 characters, so that
every q-gram also occurs exactly once as a substring of the resulting
linear sequence of length |A|^q+q-1.
'''

def de_bruijn_char_codes(alpha_size: int, q: int) -> list:
    a = [0] * alpha_size * q
    char_codes = list()

    def dbseq(t, p):
        if t > q:
            if q % p == 0:
                char_codes.extend(a[1 : p + 1])
        else:
            a[t] = a[t - p]
            dbseq(t + 1, p)
            for j in range(a[t - p] + 1, alpha_size):
                a[t] = j
                dbseq(t + 1, t)

    dbseq(1, 1)
    char_codes.extend(char_codes[:q - 1])
    return char_codes

if __name__ == '__main__':
  import sys
  import argparse
  def parse_command_line(argv):
    p = argparse.ArgumentParser(description=('generate de Bruijn sequence with '
                                             'respect to q and given alphabet'))
    p.add_argument('-q', '--qgram_length', type=int, default=3, metavar='<int>',
                    help='specify parameter q')
    alphabet_group = p.add_mutually_exclusive_group(required=True)
    p.add_argument('-t', '--transform', action='store_true', default=False,
                    help=('use integer code from 0 to alpha_size-1 for each '
                          'character'))
    alphabet_group.add_argument('-a', '--alphabet', type=str, 
                                metavar='<characters>',
                                help='specify alphabet')
    alphabet_group.add_argument('--dna', action='store_true', default=False,
                                help='use DNA alphabet')
    alphabet_group.add_argument('--protein', action='store_true', default=False,
                                help='use Protein alphabet')
    p.add_argument('--fasta', action='store_true', default=False,
                   help='output sequences in fasta format')
    return p.parse_args(argv)

  args = parse_command_line(sys.argv[1:])
  if args.dna:
    args.alphabet = 'ACGT'
  elif args.protein:
    args.alphabet = 'ADEFGHIKLMNPQRSTVWY'
  if args.transform:
    args.alphabet = ''.join(map(chr, range(len(args.alphabet))))
  char_codes = de_bruijn_char_codes(len(args.alphabet), args.qgram_length)
  seq = ''.join(args.alphabet[i] for i in char_codes)
  if args.fasta:
    print(f'>de Bruijn sequence for alphabet={args.alphabet}, '
          f'q={args.qgram_length}')
    width = 60
    for start in range(0, len(seq), width):
      print(seq[start:start + width])
  else:
    print(seq, end='')

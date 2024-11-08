#!/usr/bin/env python3

import argparse
import subprocess
import sys

def main():
    p = argparse.ArgumentParser(description='verify that unwords yields q+1 for a de Bruijn sequence')
    p.add_argument('-q', help='q-value of the de Bruijn sequence (we expect unwords of size q+1)')
    p.add_argument('inputfile', help='input file containing our sequence in FASTA format')
    args = p.parse_args()
    q = int(args.q)
    inputfile = args.inputfile
    de_bruijn_result = subprocess.getoutput(f"./unwords_mn.x --ignore_rc {inputfile}")
    if f"# length of qgrams:\t{q+1}" in de_bruijn_result:
        exit(0)
    print(de_bruijn_result, file=sys.stderr)
    raise ValueError(f"Length of unwords on de Bruijn sequence for q={q} is not {q+1}")


if __name__ == "__main__":
    main()

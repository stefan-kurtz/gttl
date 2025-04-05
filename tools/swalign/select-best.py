#!/usr/bin/env python3

import sys, re, argparse
from SEmatch import SEmatch

def show_best_hits(best,sort_by_score_only,fields_list,current_hits):
  count_best = 0
  for hit in sorted(current_hits,
                    key=lambda mdict: (mdict['score'],\
                                       0 if sort_by_score_only \
                                       else mdict['s_len'] + mdict['q_len']),\
                    reverse=True):
    l = list()
    for f in fields_list:
      l.append(hit[f])
    print('\t'.join(list(map(str,l))))
    count_best += 1
    if count_best >= best:
      return

def parse_command_line(argv):
  p = argparse.ArgumentParser(description='select given number best matches')
  p.add_argument('-a','--from_all_matches',action='store_true',default=False,
                  help=('select best matches from all matches, not from '
                        'each pair the'))
  p.add_argument('-s','--sort_by_score_only',action='store_true',default=False,
                  help=('only sort by score, but not by aligned length if '
                        'score is equal'))
  p.add_argument('best',type=int,
                  help='specify number of best hits')
  p.add_argument('matchfile',type=str,
                 help='specify input file, - means stdin')
  return p.parse_args(argv)

def main(best,from_all_matches,sort_by_score_only,matchfile):
  matchiterator = SEmatch(matchfile)
  current_hits = list()
  current_subject = None
  fields_list = None
  def convert_field(f):
    return re.sub('bit_score','bit score',re.sub(r'([sq])_',r'\1. ',f))
  for mdict in matchiterator.each():
    s_seqnum = mdict["s_seqnum"]
    if fields_list is None:
      fields_list = matchiterator.fields_get()
      assert len(fields_list) > 0
      print('# Fields: {}'.format(', '.join(map(convert_field,fields_list))))
    if not from_all_matches:
      if (current_subject is None) or (current_subject != s_seqnum):
        if len(current_hits) > 0:
          show_best_hits(best,sort_by_score_only,fields_list,current_hits)
          current_hits.clear()
        current_subject = s_seqnum
    current_hits.append(mdict)
  if len(current_hits) > 0:
    show_best_hits(best,sort_by_score_only,fields_list,current_hits)

args = parse_command_line(sys.argv[1:])
main(args.best,args.from_all_matches,args.sort_by_score_only,args.matchfile)

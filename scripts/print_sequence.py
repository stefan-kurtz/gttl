import sys

#lst{printSequence}
def print_sequence(seq,linelength,stream = sys.stdout):
  for startpos in range(0,len(seq),linelength):
    stream.write('{}\n'.format(seq[startpos:startpos+linelength]))
#lstend#

import sys

def coords_contained_in(start0,end0,start1,end1):
  if start1 <= start0 and end0 <= end1:
    return True
  return False

def coords_overlap_size(start0,end0,start1,end1):
  if start0 > start1:
    sys.stderr.write('{}: start0={} > {}=start1 not expected\n'
                      .format(sys.argv[0],start0,start1))
    exit(1)
  if end0 < start1:
    return 0
  if end0 < end1:
    return end0 - start1 + 1
  return end1 - start1 + 1

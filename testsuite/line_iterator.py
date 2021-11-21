#!/usr/bin/env python3

import sys, re

print('#!/bin/sh')
print('set -e -x')
print('TMPFILE_src=`mktemp TMP.XXXXXX` || exit 1')
print('TMPFILE_target=`mktemp TMP.XXXXXX` || exit 1')
for filename in sys.argv[1:]:
  print('cp {} ${{TMPFILE_src}}'.format(filename))
  print('./line_iterator.x ${TMPFILE_src} | diff - ${TMPFILE_src}')
  for copy_num in range(2,4):
    file_list = ' '.join(['${TMPFILE_src}'] * copy_num)
    print('cat {} > ${{TMPFILE_target}}'.format(file_list))
    print('./line_iterator.x {} | diff - ${{TMPFILE_target}}'.format(file_list))

print('rm -f ${TMPFILE_src} ${TMPFILE_target}')

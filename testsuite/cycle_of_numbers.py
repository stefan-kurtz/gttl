#!/usr/bin/env python3

excluded = {ord(cc) for cc in 'ACGTNacgtn'}

for i in range(255+1):
  if i not in excluded:
    print(i)

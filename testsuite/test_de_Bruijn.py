#!/usr/bin/env python3
import sys
sys.path.append(r'../scripts')
from de_Bruijn import de_bruijn_char_codes

for alpha_size in [2, 3, 4, 20]:
  for q in range(1, 11):
    if alpha_size ** q > 2_000_000:
      break
    seq = de_bruijn_char_codes(alpha_size, q)
    assert len(seq) == alpha_size ** q + q - 1
    qgrams = {tuple(seq[i:i + q]) for i in range(len(seq) - q + 1)}
    assert len(qgrams) == alpha_size ** q, (alpha_size, q)
    print(f'de_Bruijn_sequence for alpha_size={alpha_size}, q={q}: ok')

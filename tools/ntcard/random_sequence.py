from random import random

class RandomSequence:
  def __init__(self,characters,probs):
    assert len(characters) > 0
    assert len(characters) == len(probs)
    assert sum(probs) == 1.0
    self.cummulative = list()
    psum = 0
    for p in probs:
      psum += p
      self.cummulative.append(psum)
    assert self.cummulative[-1] == 1.0
    assert len(self.cummulative) == len(probs)
    self.characters = characters
  def alphasize(self):
    return len(self.characters)
  def rchar(self):
    r = random()
    assert r >= 0 and r <= 1.0
    for idx, c in enumerate(self.cummulative):
      if r < c:
        return self.characters[idx]
    assert False
    return None
  def rseq(self,this_length):
    return ''.join([self.rchar() for _ in range(this_length)])

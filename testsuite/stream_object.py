import sys, re, gzip

class StreamObject:
  def __init__(self,input_object):
    self.decode = False
    self.do_close = False
    if isinstance(input_object,str):
      if input_object == '-':
        self.stream = sys.stdin
        self.matchinput_message = 'STDIN'
      else:
        self.matchinput_message = input_object
        try:
          if re.search(r'\.gz$',input_object):
            self.stream = gzip.open(input_object,'r')
            self.decode = True
          else:
            self.stream = open(input_object,'r')
            self.do_close = True
        except IOError as err:
          sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
          exit(1)
    else:
      self.stream = input_object
      self.matchinput_message = 'GENERATOR'
  def __iter__(self):
    for line in self.stream:
      if self.decode:
        yield line.decode('utf-8').rstrip()
      else:
        yield line.rstrip()
  def __del__(self):
    if self.do_close:
      self.stream.close()

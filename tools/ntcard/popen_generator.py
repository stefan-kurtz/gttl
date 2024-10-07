import sys, shlex, subprocess

'''the following function executes the given command given as string
   cmd_line. Each line in the output is captures and delivered by
   a yield statement. So the typical use is to execute the function
   in a for-loop like this:
   for line in popen_generator('ls -l'):
     values = line.rstrip().split()
     ...
'''

def popen_generator(cmd_line,with_stdin=True, with_stderr=False):
  cmd_args = shlex.split(cmd_line)
  thispipe = subprocess.Popen(cmd_args,
                              stdin=subprocess.PIPE if with_stdin else None,
                              stderr=subprocess.PIPE if with_stderr else None,
                              stdout=subprocess.PIPE)
  if with_stderr:
    out = thispipe.stderr.read()
  else:
    # fill pipe to program delivering input using converted_input.encode()
    out, _ = thispipe.communicate()
  exit_code = thispipe.wait()
  if exit_code != 0:
    sys.stderr.write('{}: FAILURE: {}\n'.format(sys.argv[0],cmd_line))
    exit(1)
  line = list()
  for cc in out.decode():
    if cc == '\n':
      yield ''.join(line)
      line.clear()
    else:
      line.append(cc)

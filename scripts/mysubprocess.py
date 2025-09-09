import sys, re, subprocess, shlex

def stream2string(stream):
  s = str()
  for cc in stream.decode():
    s += cc
  return s.rstrip()

def remove_suppression(s):
  return re.sub(r'\n----------.*$','',s,flags=re.S)

def mysubprocess(cmd_line):
  cmd_args = shlex.split(cmd_line)
  thispipe = subprocess.Popen(cmd_args,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
  stdout_encoded_stream, stderr_encoded_stream = thispipe.communicate(0)
  rc = thispipe.returncode
  stdout_str = stream2string(stdout_encoded_stream)
  stderr_err = stream2string(stderr_encoded_stream)
  return rc, stdout_str, remove_suppression(stderr_err)

def normalize_err_msg(s):
  if isinstance(s, str):
    return re.sub(r'^[^:]+: ','',s)
  return s

def mysubprocess_expect(cmd_line,expected_err_code,expected_err_msg = None):
  rc, stdout_str, stderr_str = mysubprocess(cmd_line)
  if rc != expected_err_code:
    sys.stderr.write((f'{sys.argv[0]}: cmd_line="{cmd_line}", err_code = {rc} '
                      f'!= {expected_err_code} = expected_err_code\n'))
    exit(1)
  if expected_err_msg:
    n_stderr_str = normalize_err_msg(stderr_str)
    n_expected_err_msg = normalize_err_msg(expected_err_msg)
    if isinstance(expected_err_msg,set):
      for word in expected_err_msg:
        if word not in n_stderr_str:
          sys.stderr.write((f'err_code={rc}: cmd_line="{cmd_line}",\n'
                            f'stderr_str = "{n_stderr_str}" '
                            f'does not contain the keyword {word}\n'))
          exit(1)
    elif not n_stderr_str.startswith(n_expected_err_msg):
      sys.stderr.write((f'{sys.argv[0]}: cmd_line="{cmd_line}",\n'
                        f'stderr_str =\n"{n_stderr_str}" !=\n'
                        f'"{n_expected_err_msg}"\n= expected_err_msg\n'))
      if len(n_stderr_str) != len(n_expected_err_msg):
        sys.stderr.write((f'diagnosis: len(n_stderr_str) = {len(n_stderr_str)} '
                          f'!= {len(n_expected_err_msg)} = '
                          'len(n_expected_err_msg)\n'))
      else:
        diff_pos = [i for i in range(len(n_stderr_str))
                      if n_stderr_str[i] != n_expected_err_msg[i]]
        assert len(diff_pos) > 0
        first_diff_pos = diff_pos[0]
        sys.stderr.write((f'diagnosis: n_stderr_str[{first_diff_pos}] = '
                          f'{n_stderr_str[first_diff_pos]} != '
                          f'{n_expected_err_msg[first_diff_pos]} = '
                          f'n_expected_err_msg[{first_diff_pos}]\n'))
      exit(1)

#!/usr/bin/env python3

'''
Eingabe: Dateinamen mit Programmcode (z.B. in Python oder C). Finde
Formatierungsfehler in der Datei und breche beim ersten solchen Fehler
ab. Die folgendne Fehler werden mit einer Zeilennummer und einer entsprechenden
Meldung berichtet und das Programm bricht mit Fehlercode 1 ab:
- Zeilen, die mehr als 80 Zeichen lang sind

optional:
- Option -t: Vorkommen von Tabulatoren
- Option -w: white space am Ende einer Zeile
- Nur sinnvoll f"ur Python:
  Option -i: Einrueckung von Zeilen, die nicht mit einem
  der Zeichen ."'# beginnen und nicht durch 4 oder 2 ganzzahlig teilbar sind

Da Programm kann beliebig viele Dateien als Argumente verarbeiten.
Ein Argument - wird als sys.stdin interpretiert.
'''

import sys, re, argparse

def parse_command_line(argv):
  p = argparse.ArgumentParser()
  p.add_argument('-w','--trailing_white_space',action='store_true',
                 default=False,help='check for trailing white spaces')
  p.add_argument('-i','--indent_level',action='store_true',
                 default=False,help='check for indentation level')
  p.add_argument('-t','--no_tabulator',action='store_true',
                 default=False,help='check for occurrence of tabulator')
  p.add_argument('inputfile',nargs='+',
                  help='specify input file, - means stdin')
  return p.parse_args(argv)

def errormsg(error_list,filename,linenum,msg):
  error_list.append('file {}, line {}: {}'.format(filename,linenum,msg))

def code_check(error_list,trailing_white_space,filename):
  if filename == '-':
    stream = sys.stdin
  else:
    #sys.stderr.write('{}: open {}\n'.format(sys.argv[0],filename))
    try:
      stream = open(filename)
    except IOError as err:
      sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
      exit(1)
  linenum = 0
  for line in stream:
    if len(error_list) >= 10:
      break
    line = line.rstrip('\n')
    linenum += 1
    if args.no_tabulator and re.search(r'\t',line):
      errormsg(error_list,filename,linenum,'line contains tabulator symbol')
    if trailing_white_space and re.search(r'\s+$',line):
      errormsg(error_list,filename,linenum,'trailing white space')
    if len(line) > 80 and (not line.startswith('#include')):
      errormsg(error_list,filename,linenum,'line longer than 80 characters')
    if args.indent_level:
      m = re.search(r'(^\s+)(\S)',line)
      if m:
        indent = m.group(1)
        indent_size = len(indent)
        if indent_size % 2 != 0 and indent_size % 4 != 0:
          first_char = m.group(2)
          if not (first_char in '."#\''):
            errormsg(error_list,filename,linenum,\
                     'indentation level must be multiples of 2 or 4')
  stream.close()

args = parse_command_line(sys.argv[1:])

error_list = list()
for filename in args.inputfile:
  code_check(error_list,args.trailing_white_space,filename)

for err_msg in error_list:
  sys.stderr.write('{}: {}\n'.format(sys.argv[0],err_msg))

if len(error_list) > 0:
  exit(1)

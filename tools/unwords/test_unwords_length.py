#!/usr/bin/env python3
import sys, subprocess, re, os

def show_unwords_length(filename, qmax):
    output = subprocess.getoutput("./unwords_mn.x {} {}".format(qmax,filename))
    file_size = re.search(r'total size of inputfiles \(B\):\s+(\d+)', output)
    est_qgram_length = re.search(r'estimated length of unwords:\s+(\d+)',
                                 output)
    qgram_length = re.search(r'length of qgrams:\s+(\d+)', output)
    if file_size and est_qgram_length and qgram_length:
        print("{}\t{}\t{}".format(file_size.group(1), qgram_length.group(1),
                                  est_qgram_length.group(1)))
        return True
    else:
        return False

def main():
    if len(sys.argv) != 2:
        sys.stderr.write('{}: <testdata path>\n'.format(sys.argv[0]))
        exit(1)

    qmax = 16
    path = sys.argv[1]
    files = " ".join(next(os.walk(path))[2])
    filenames = ""
    for filename in re.findall(r'\S+\.fna', files):
        filename = path+"/"+filename
        if show_unwords_length(filename, qmax): # Dateien einzeln
            filenames += filename+" "
            show_unwords_length(filenames, qmax) # Dateien zsm

if __name__ == '__main__':
    main()
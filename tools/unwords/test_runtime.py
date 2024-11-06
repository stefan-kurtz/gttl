#!/usr/bin/env python3
import sys, re, argparse, subprocess

# Ausgabe Computerdaten
def print_computer_data():
    print('# Daten des verwendeten Computers:')

    output = subprocess.getoutput('lscpu')

    s = re.search(r'Model name:\s+(.+)', output)
    print('# Modellname\t{}'.format(s.group(1)))

    s = re.search(r'BogoMIPS:\s+(.+)', output)
    print('# BogoMIPS\t{}'.format(s.group(1)))

    print('# Caches:')
    for s in re.findall(r'(.+Cache:)\s+(.+)', output):
        print('# {}\t{}'.format(s[0], s[1]))

    print()

# Ausgabe Laufzeittests
def runtime_test_unwords_mn(testdata_descr, testdata,
                            qmax_descr, qmax,
                            RC_descr, RC):
    output = subprocess.getoutput('./unwords_mn.x -s {} {} {}'
                                  .format(RC, qmax, testdata))
    s = re.search(r'compute unwords by a binary search \(ms\):\s+(\d+)',output)
    if s:
        print('unwords_mn\t{}\t{}\t{}\t{} ms'
              .format(testdata_descr, qmax_descr, RC_descr, s.group(1)))

def runtime_test_unwords_bits(testdata_descr, testdata):
    output = subprocess.getoutput('time ../../unwords-bits/Unwords-bits-src/'
                                  'unwords-bits {}'
                                  .format(testdata))
    s = re.search(r'(\d+):(\d+).(\d+)elapsed', output)
    if s:
        runtime = int(s.group(1))*60000+int(s.group(2))*1000+int(s.group(3))*10
        # Minuten*60000 + Sekunden*1000 + (Millisekunden/10)*10
        print('unwords-bits\t{}\t16\tmit RC\t{} ms'
              .format(testdata_descr, runtime))

def parse_command_line(args):
    p = argparse.ArgumentParser(description='generate runtime tests')
    p.add_argument('--computer_data',action='store_true',default=False,
                   help='print computer data')
    p.add_argument('--unwords_mn',action='store_true',default=False,
                   help='test runtime of unwords_mn')
    p.add_argument('--unwords_bits',action='store_true',default=False,
                   help='test runtime of unwords-bits')
    p.add_argument('--protein',action='store_true',default=False,
                   help='use protein testdata')
    return p.parse_args(args)

def main():
    args = parse_command_line(sys.argv[1:])
    testdata_list  = [['at1MB', '../../gttl/testdata/at1MB.fna'],\
                      ['vaccg', '../../gttl/testdata/vaccg.fna'],\
                      ['ychrIII', '../../gttl/testdata/ychrIII.fna'],\
                      ['MusI', '../../testdata/Mus_musculus.GRCm39.dna.'
                               'chromosome.1.fa'],\
                      ['human', '../../testdata/GCA_009914755.4_T2T-CHM13v2.0'
                                '_genomic.fna']]
    qmax_list = [10, 16]
    if args.protein:
        testdata_list = [['human1', '../../testdata/human.1.protein.faa'],
                         ['human2', '../../testdata/human.2.protein.faa'],
                         ['human3', '../../testdata/human.3.protein.faa'],
                         ['human', '../../testdata/human.1.protein.faa '\
                                   '../../testdata/human.2.protein.faa '\
                                   '../../testdata/human.3.protein.faa'],
                         ['mouse1', '../../testdata/mouse.1.protein.faa'],
                         ['mouse2', '../../testdata/mouse.2.protein.faa'],
                         ['mouse', '../../testdata/mouse.1.protein.faa '\
                                   '../../testdata/mouse.2.protein.faa'],
                         ['frog1', '../../testdata/frog.1.protein.faa'],
                         ['pig1', '../../testdata/pig.1.protein.faa'],
                         ['rat1', '../../testdata/rat.1.protein.faa'],
                         ['rat2', '../../testdata/rat.2.protein.faa'],
                         ['rat', '../../testdata/rat.1.protein.faa '\
                                 '../../testdata/rat.2.protein.faa']]

    if args.computer_data:
        print_computer_data()

    if args.unwords_mn or args.unwords_bits:
        print('# Programm\tDaten\tqmax\tRC\tLaufzeit')

        if args.unwords_mn:
            for testdata in testdata_list:
                for qmax in qmax_list:
                    runtime_test_unwords_mn(testdata[0], testdata[1],
                                            qmax, '-q {} -i'.format(qmax),
                                            'ohne RC', '')
                    runtime_test_unwords_mn(testdata[0], testdata[1],
                                            qmax, '-q {}'.format(qmax),
                                            'mit RC', '-i')
                runtime_test_unwords_mn(testdata[0], testdata[1],
                                        'est.', '',
                                        'ohne RC', '-i')
                runtime_test_unwords_mn(testdata[0], testdata[1],
                                        'est.', '',
                                        'mit RC', '')

        if args.unwords_bits:
            for testdata in testdata_list:
                    runtime_test_unwords_bits(testdata[0], testdata[1])

if __name__ == '__main__':
    main()

import sys
import gzip
import bz2
import os
import time
import subprocess
import math
import multiprocessing as mp
from optparse import OptionParser

import numpy as np

from utils import *
import AAF

def main():
    parser = OptionParser()
    parser.add_option("-k", dest="k_len", type=int, default=25,
                      help="k for reconstruction, default = 25")
    parser.add_option("--k_s", dest="k_s_len", type=int, default=25,
                      help="k for reads selection, default = 25")
    parser.add_option("-n", dest="filter", type=int, default=1,
                      help="k-mer filtering threshold, default = 1")
    parser.add_option("-d", dest="data_dir", default='data',
                      help="directory containing  data, default = data/")
    parser.add_option("-G", dest="mem_size", type=int, default=4,
                      help="total memory limit (GB), default = 4")
    parser.add_option("-t", dest="n_threads", type=int, default=1,
                      help="number of threads to use, default = 1")
    parser.add_option("-l", dest="long", action='store_true',
                      help="use fitch_kmerX_long instead of fitch_kmerX")

    options, args = parser.parse_args()

    n = options.filter
    mem_size = options.mem_size
    n_threads = options.n_threads
    k_l = options.k_len
    k_s = options.k_s_len
    mem_per_thread = int(options.mem_size / float(n_threads))
    data_dir = options.data_dir

    if not mem_per_thread:
        print('Not enough memory, decrease n_threads or increase mem_size!')
        sys.exit()

    if not os.path.isdir(data_dir):
        print('Cannot find data directory {}!'.format(data_dir))
        sys.exit(2)


    if k_l > 25:
        if os.system('which kmer_countx > /dev/null'):
            kmer_count = './kmer_countx'
            if not is_exe(kmer_count):
                print('kmer_countx not found!')
                sys.exit(1)
        else:
            kmer_count = 'kmer_countx'

    else:
        if os.system('which kmer_count > /dev/null'):
            kmer_count = './kmer_count'
            if not is_exe(kmer_count):
                print('kmer_count not found!')
                sys.exit(1)
        else:
            kmer_count = 'kmer_count'

    if os.system('which kmer_merge > /dev/null'):
        filt = './kmer_merge'
        if not is_exe(filt):
            print('kmer_merge not found!')
            sys.exit(1)
    else:
        filt = 'kmer_merge'

    if os.system('which ReadsSelector > /dev/null'):
        ReadsSelector = './ReadsSelector'
        if not is_exe(filt):
            print('ReadsSelector not found!')
            sys.exit(1)
    else:
        ReadsSelector = 'ReadsSelector'

    if os.system('which fitch_kmerX > /dev/null'):
        if options.long:
            fitch = './fitch_kmerX_long'
        else:
            fitch = './fitch_kmerX'
        if not is_exe(fitch):
            print(fitch + ' not found!')
            sys.exit()
    else:
        if options.long:
            fitch = 'fitch_kmerX_long'
        else:
            fitch = 'fitch_kmerX'


    selection_dir = '{}_k_s{}_pairwise'.format(
        os.path.basename(data_dir.rstrip('/')), k_s)

    if os.path.exists('./' + selection_dir):
        command = 'rm -r {}'.format(selection_dir)
        os.system(command)
    command = 'mkdir {}'.format(selection_dir)
    os.system(command)



    samples = AAF.aaf_kmer_count(data_dir, k_s, n, n_threads, mem_per_thread)

    # build distance matrix
    sn = len(samples)
    dist = [[0] * sn for i in range(sn)]
    for i in range(sn):
        for j in range(i + 1, sn):
            command = []
            command.append(
                'mkdir {}/{}_{}'.format(selection_dir, samples[i], samples[j]))
            command.append(
                '{} -k s -c -d 0 -A A -B A {}.pkdat.gz {}.pkdat.gz | cut -f 1 > test.kmer'.format(filt, samples[i], samples[j]))
            command.append('{} -k test.kmer -fa 1 -o {}/{}_{}/{} -s {}/{}/*'
                           .format(ReadsSelector, selection_dir, samples[i], samples[j],
                                   samples[i], data_dir, samples[i]))
            command.append('{} -k test.kmer -fa 1 -o {}/{}_{}/{} -s {}/{}/*'
                           .format(ReadsSelector, selection_dir, samples[i], samples[j],
                                   samples[j], data_dir, samples[j]))
            for comm in command:
                print(comm)
                os.system(comm)

            n_total = []
            command = '{} -l {} -n {} -G {} -o {}_temp.pkdat -f FA -i {}/{}_{}/{}.*'.format(
                kmer_count, k_l, n, mem_size / 2, samples[i], selection_dir, samples[i], samples[j], samples[i])
            output = subprocess.check_output(
                command, shell=True, stderr=subprocess.STDOUT)
            n_total.append(float(output.decode('ascii').split()[1]))
            command = '{} -l {} -n {} -G {} -o {}_temp.pkdat -f FA -i {}/{}_{}/{}.*'.format(
                kmer_count, k_l, n, mem_size / 2, samples[j], selection_dir, samples[i], samples[j], samples[j])
            output = subprocess.check_output(
                command, shell=True, stderr=subprocess.STDOUT)
            n_total.append(float(output.decode('ascii').split()[1]))
            command = "{} -k s -c -d '0' -a 'T,M,F' {}_temp.pkdat {}_temp.pkdat | wc -l".format(
                filt, samples[i], samples[j])
            output = subprocess.check_output(
                command, shell=True, stderr=subprocess.STDOUT)
            n_shared = int(output.decode('ascii').split()[0])
            if n_shared == 0:
                distance = 1
            else:
                distance = (-1.0 / k_l) * math.log(n_shared / min(n_total))
            dist[j][i] = dist[i][j] = distance

    os.system('rm *.pkdat*')

    # construct phylogeny
    try:
        in_file = open('in_file', 'w')
    except IOError:
        print('Cannot open in_file for writing')
        sys.exit()

    in_file.write('{} {}'.format(sn, sn))
    names = {}
    for i in range(sn):
        lsl = len(samples[i])
        if lsl >= 10:
            ssl = samples[i][:10]
            appendix = 1
            while ssl in names:
                if appendix < 10:
                    ssl = samples[i][:9] + str(appendix)
                elif appendix > 9:
                    ssl = samples[i][:8] + str(appendix)
                appendix += 1
        else:
            ssl = samples[i] + ' ' * (10 - lsl)
        names[ssl] = samples[i]
        in_file.write('\n{}'.format(ssl))
        for j in range(sn):
            in_file.write('\t{}'.format(dist[i][j]))

    in_file.close()

    # run fitch_kmer
    print('{} building tree'.format(time.strftime("%c")))
    if os.path.exists("./outfile"):
        os.system("rm -f outfile outtree")
    command = 'printf "K\n{}\nY" | {} > /dev/null'.format(int(k_l), fitch)
    os.system(command)
    fh0 = open('outtree', 'rt')
    fh1 = open(selection_dir + '.tre', 'wt')

    for line in fh0:
        for key in names:
            key_new = key.rstrip() + ":"
            if key_new in line:
                newline = line.replace(key_new, names[key].rstrip() + ":", 1)
                line = newline
        fh1.write(line)
    fh0.close()
    fh1.close()
    command = 'mv in_file {}.dist'.format(selection_dir)
    os.system(command)

    os.system('rm -f outfile outtree')

    print('{} end'.format(time.strftime("%c")))


if __name__ == '__main__':
    main()

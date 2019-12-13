import sys
import os
import time
import math
import psutil
from concurrent.futures import ProcessPoolExecutor as PPE
import numpy as np

from utils import *


def aaf_kmercount(data_dir, k, n, n_threads, mem_per_thread):
    if k > 25:
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

    samples = []
    for filename in os.listdir(data_dir):
        if os.path.isdir(os.path.join(data_dir, filename)):
            samples.append(filename)
        else:
            if not filename.startswith('.'):
                sample = filename.split(".")[0]
                if sample in samples:
                    sample = filename.split(".")[0] + filename.split(".")[1]
                    if sample in samples:
                        print('Error! Redundant sample or file names! Aborting!')
                        sys.exit(3)
                os.system("mkdir {}/{}".format(data_dir, sample))
                os.system("mv {}/{} {}/{}/".format(data_dir,
                                                   filename, data_dir, sample))
                samples.append(sample)
    samples.sort()
    print(time.strftime('%c'))
    print('SPECIES LIST:')
    for sample in samples:
        print(sample)

    jobList = []
    for sample in samples:
        outFile = '{}.pkdat.gz'.format(sample)
        command = '{} -l {} -n {} -G {} -o {} -f '.format(kmer_count, k, n,
                                                          mem_per_thread, outFile)
        c = ''

        for input_file in os.listdir(os.path.join(data_dir, sample)):
            input_file = os.path.join(data_dir, sample, input_file)
            handle = smart_open(input_file)
            firstChar = handle.read(1)
            if firstChar == '@':
                seqFormat = 'FQ'
            elif firstChar == '>':
                seqFormat = 'FA'
            else:
                print(
                    'Error, file {} is not FA or FQ format. Aborting!'.format(input_file))
                sys.exit(3)
            c += " -i '{}'".format(input_file)
        command += '{}{}> {}.wc'.format(seqFormat, c, sample)
        jobList.append(command)
    # Run jobs
    with PPE(max_workers=n_threads) as executor:
        executor.map(run_command, jobList)
    return samples


def aaf_dist(datfile, count_file, n_threads, samples, kl, long=False):
    if os.system('which fitch_kmerX > /dev/null'):
        if long:
            fitch = './fitch_kmerX_long'
        else:
            fitch = './fitch_kmerX'
        if not is_exe(fitch):
            print(fitch + ' not found. Make sure it is in your PATH or the')
            print('current directory, and that it is executable')
            sys.exit()
    else:
        if long:
            fitch = 'fitch_kmerX_long'
        else:
            fitch = 'fitch_kmerX'
    # process the .dat.gz file
    try:
        iptf = smart_open(datfile, 'rt')
    except IOError:
        print('Cannot open file', data_file)
        sys.exit()

    if not os.path.isfile(count_file):
        print('Cannot find file', count_file)
        sys.exit()

    try:
        total = open(count_file, 'rt')
    except IOError:
        print('Cannot open file', count_file)
        sys.exit()

    try:
        in_file = open('in_file', 'wt')
    except IOError:
        print('Cannot open in_file for writing')
        sys.exit()

    # Read header
    sl = []  # species list
    line = iptf.readline()
    ll = line.split()
    if kl != float(ll[1]):
        # kmer length
        print("k in kmer table file not same as k given to aaf_dist!")
        sys.exit()
    while True:
        line = iptf.readline()
        if line.startswith('#-'):
            continue
        elif line.startswith('#sample'):
            ll = line.split()
            sl.append(ll[1])
        else:
            break
    if sl != samples:
        print("sample list in  kmer table file not same as k given to aaf_dist!")
        sys.exit()
    s_n = len(samples)  # species number
    n_share = [[0] * s_n for i in range(s_n)]

    cpu_num = psutil.cpu_count()
    line = iptf.readline()
    line_size = sys.getsizeof(line)
    chunk_length = int(1024 ** 3 / cpu_num / line_size)
    print('chunk_length = {}'.format(chunk_length))
    while True:
        lines = []
        for nLines in range(chunk_length):
            if not line:
                break
            lines.append(line)
            line = iptf.readline()
        if not lines:
            break
        with PPE(max_workers=cpu_num) as executor:
            for result in executor.map(count_shared_single, lines):
                for i in range(s_n):
                    for j in range(i + 1, s_n):
                        n_share[i][j] += result[i][j]

    iptf.close()

    # Compute distance matrix
    n_total = [0.0] * s_n

    for i in range(s_n):
        n_total[i] = float(total.readline().split()[1])
    dist = [[0] * s_n for i in range(s_n)]

    for i in range(s_n):
        for j in range(i + 1, s_n):
            min_total = min(n_total[i], n_total[j])
            if n_share[i][j] == 0:
                dist[j][i] = dist[i][j] = 1
            else:
                distance = (-1 / float(kl) * math.log(n_share[i][j] / min_total))
                dist[j][i] = dist[i][j] = distance
                n_share[j][i] = n_share[i][j]

    total.close()

    # Write in_file
    in_file.write('{} {}'.format(s_n, s_n))
    names = {}
    for i in range(s_n):
        lsl = len(sl[i])
        if lsl >= 10:
            ssl = sl[i][:10]
        appendix = 1
        while ssl in names:
            ssl = ssl[:-len(str(appendix))] + str(appendix)
            appendix += 1
        if lsl < 10:
            ssl = sl[i] + ' ' * (10 - lsl)
        names[ssl] = sl[i]
        in_file.write('\n{}'.format(ssl))
        for j in range(s_n):
            in_file.write('\t{}'.format(dist[i][j]))

    in_file.close()

    # Run fitch_kmer
    print('{} building tree'.format(time.strftime("%c")))
    if os.path.exists("./outfile"):
        os.system("rm -f outfile outtree")
    command = 'printf "K\n{}\nY" | {} > /dev/null'.format(int(kl), fitch)
    os.system(command)
    fh = open('outtree', 'rt')
    fh1 = open(datfile.split('.')[0] + '.tre', 'wt')

    for line in fh:
        for key in names:
            key_new = key.rstrip() + ":"
            if key_new in line:
                newline = line.replace(key_new, names[key].rstrip() + ":", 1)
                line = newline
        fh1.write(line)
    fh.close()
    fh1.close()
    command = 'mv in_file {}.dist'.format(datfile.split('.')[0])
    os.system(command)

    os.system('rm -f outfile outtree')

    print('{} end'.format(time.strftime("%c")))

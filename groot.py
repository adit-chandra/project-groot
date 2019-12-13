import sys
import gzip
import bz2
import os
import time
import subprocess
import math
import multiprocessing as mp
from optparse import OptionParser

from utils import *
import AAF


def main():
    parser = OptionParser()
    parser.add_option("-k", dest="kLen", type=int, default=25,
                      help="k for reconstruction, default = 25")
    parser.add_option("--ks", dest="ksLen", type=int, default=25,
                      help="k for reads selection, default = 25")
    parser.add_option("-n", dest="filter", type=int, default=1,
                      help="k-mer filtering threshold, default = 1")
    parser.add_option("-d", dest="data_dir", default='data',
                      help="directory containing the data, default = data/")
    parser.add_option("--mem", dest="mem_size", type=int, default=4,
                      help="total memory limit (in GB), default = 4")
    parser.add_option("-t", dest="n_threads", type=int, default=1,
                      help="number of threads to use, default = 1")
    parser.add_option("-l", dest="long", action='store_true',
                      help="use fitch_kmerX_long instead of fitch_kmerX")

    (options, args) = parser.parse_args()

    n = options.filter
    mem_size = options.mem_size
    n_threads = options.n_threads
    kl = options.kLen
    ks = options.ksLen
    mem_per_thread = int(options.mem_size / float(n_threads))
    data_dir = options.data_dir

    if not mem_per_thread:
        print('Not enough memory, decrease n_threads or increase mem_size')
        sys.exit()

    if not os.path.isdir(data_dir):
        print('Cannot find data directory {}'.format(data_dir))
        sys.exit(2)

    # kmer_count
    if kl > 25:
        if os.system('which kmer_countx > /dev/null'):
            kmer_count_bin = './kmer_countx'
            if not is_exe(kmer_count_bin):
                print('kmer_countx not found!')
                sys.exit(1)
        else:
            kmer_count_bin = 'kmer_countx'

    else:
        if os.system('which kmer_count > /dev/null'):
            kmer_count_bin = './kmer_count'
            if not is_exe(kmer_count_bin):
                print('kmer_count not found!')
                sys.exit(1)
        else:
            kmer_count_bin = 'kmer_count'

    # kmer_merge
    if os.system('which kmer_merge > /dev/null'):
        filt_bin = './kmer_merge'
        if not is_exe(filt_bin):
            print('kmer_merge not found!')
            sys.exit(1)
    else:
        filt_bin = 'kmer_merge'

    # ReadsSelector
    if os.system('which ReadsSelector > /dev/null'):
        reads_select_bin = './ReadsSelector'
        if not is_exe(filt_bin):
            print('ReadsSelector not found!')
            sys.exit(1)
    else:
        reads_select_bin = 'ReadsSelector'

    # fitch
    if os.system('which fitch_kmerX > /dev/null'):
        if options.long:
            filtch_bin = './fitch_kmerX_long'
        else:
            filtch_bin = './fitch_kmerX'
        if not is_exe(filtch_bin):
            print(filtch_bin + ' not found')
            sys.exit()
    else:
        if options.long:
            filtch_bin = 'fitch_kmerX_long'
        else:
            filtch_bin = 'fitch_kmerX'


    selection_dir = '{}_ks{}_pairwise'.format(
        os.path.basename(data_dir.rstrip('/')), ks)

    if os.path.exists('./' + selection_dir):
        command = 'rm -r {}'.format(selection_dir)
        os.system(command)
    command = 'mkdir {}'.format(selection_dir)
    os.system(command)



    samples = AAF.aaf_kmer_count(data_dir, ks, n, n_threads, mem_per_thread)

    # build distance matrix
    sn = len(samples)
    dist = [[0] * sn for i in range(sn)]
    for i in range(sn):
        for j in range(i + 1, sn):
            command = []
            command.append(
                'mkdir {}/{}_{}'.format(selection_dir, samples[i], samples[j]))
            command.append(
                '{} -k s -c -d 0 -A A -B A {}.pkdat.gz {}.pkdat.gz | cut -f 1 > test.kmer'.format(filt_bin, samples[i], samples[j]))
            command.append('{} -k test.kmer -fa 1 -o {}/{}_{}/{} -s {}/{}/*'
                           .format(reads_select_bin, selection_dir, samples[i], samples[j],
                                   samples[i], data_dir, samples[i]))
            command.append('{} -k test.kmer -fa 1 -o {}/{}_{}/{} -s {}/{}/*'
                           .format(reads_select_bin, selection_dir, samples[i], samples[j],
                                   samples[j], data_dir, samples[j]))
            for comm in command:
                print(comm)
                os.system(comm)

            # kmer_count
            ntotal = []
            command = '{} -l {} -n {} -G {} -o {}_temp.pkdat -f FA -i {}/{}_{}/{}.*'.format(
                kmer_count_bin, kl, n, mem_size / 2, samples[i], selection_dir, samples[i], samples[j], samples[i])
            output = subprocess.check_output(
                command, shell=True, stderr=subprocess.STDOUT)
            ntotal.append(float(output.decode('ascii').split()[1]))
            command = '{} -l {} -n {} -G {} -o {}_temp.pkdat -f FA -i {}/{}_{}/{}.*'.format(
                kmer_count_bin, kl, n, mem_size / 2, samples[j], selection_dir, samples[i], samples[j], samples[j])
            output = subprocess.check_output(
                command, shell=True, stderr=subprocess.STDOUT)
            ntotal.append(float(output.decode('ascii').split()[1]))

            # kmer_merge
            command = "{} -k s -c -d '0' -a 'T,M,F' {}_temp.pkdat {}_temp.pkdat | wc -l".format(
                filt_bin, samples[i], samples[j])
            output = subprocess.check_output(
                command, shell=True, stderr=subprocess.STDOUT)
            nshared = int(output.decode('ascii').split()[0])
            if nshared == 0:
                distance = 1
            else:
                distance = (-1.0 / kl) * math.log(nshared / min(ntotal))
            dist[j][i] = dist[i][j] = distance

    os.system('rm *.pkdat*')

    # construct the tree
    try:
        infile = open('infile', 'w')
    except IOError:
        print('Cannot open infile for writing')
        sys.exit()

    infile.write('{} {}'.format(sn, sn))
    namedic = {}
    for i in range(sn):
        lsl = len(samples[i])
        if lsl >= 10:
            ssl = samples[i][:10]
            appendix = 1
            while ssl in namedic:
                if appendix < 10:
                    ssl = samples[i][:9] + str(appendix)
                elif appendix > 9:
                    ssl = samples[i][:8] + str(appendix)
                appendix += 1
        else:
            ssl = samples[i] + ' ' * (10 - lsl)
        namedic[ssl] = samples[i]
        infile.write('\n{}'.format(ssl))
        for j in range(sn):
            infile.write('\t{}'.format(dist[i][j]))

    infile.close()

    # fitch_kmer
    print('{} building tree'.format(time.strftime("%c")))
    if os.path.exists("./outfile"):
        os.system("rm -f outfile outtree")
    command = 'printf "K\n{}\nY" | {} > /dev/null'.format(int(kl), filtch_bin)
    os.system(command)
    fh1 = open('outtree', 'rt')
    fh2 = open(selection_dir + '.tre', 'wt')

    for line in fh1:
        for key in namedic:
            key_new = key.rstrip() + ":"
            if key_new in line:
                newline = line.replace(key_new, namedic[key].rstrip() + ":", 1)
                line = newline
        fh2.write(line)
    fh1.close()
    fh2.close()
    command = 'mv infile {}.dist'.format(selection_dir)
    os.system(command)

    os.system('rm -f outfile outtree')

    print('{} end'.format(time.strftime("%c")))

if __name__ == '__main__':
    main()

import sys
import os
import time
import math
import psutil
from concurrent.futures import ProcessPoolExecutor as PPE

import numpy as np


def run_command(command):
    print(command)
    print(time.strftime('%c'))
    try:
        os.system(command)
    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise


def smart_open(filename, mode='rt'):
    import gzip
    import bz2
    if filename.endswith('gz'):
        return gzip.open(filename, mode)
    elif filename.endswith('bz2'):
        return bz2.BZ2File(filename, mode)
    else:
        return open(filename, mode)


def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def present(x, n=1):
    if int(x) >= n:
        return '1'
    else:
        return '0'


def count_total(lines):
    line_list = []
    for line in lines:
        line_list.append([int(present(i)) for i in line.split()[1:]])
    line_total = np.sum(line_list, axis=0)
    return line_total


def count_shared_single(line):  # count nshare only, for shared kmer table
    line = line.split()
    if line[0][0].isdigit():
        sn = len(line)
        flag = 'd'
    else:
        sn = len(line) - 1
        flag = 'k'
    shared = [[0] * sn for i in range(sn)]
    if flag == 'k':
        line = line[1:]
    line = [int(i) for i in line]
    for i in range(sn):
        for j in range(i + 1, sn):
            if line[i] * line[j] != 0:
                shared[i][j] += 1
    return shared


def count_total_shared(lines, sn):
    line_list = []
    shared = [[0] * sn for i in range(sn)]
    for line in lines:
        line = line.split()
        if len(line) == sn + 1:
            line = line[1:]
        line = [int(i) for i in line]
        line_list.append([int(present(i)) for i in line])
        for i in range(sn):
            for j in range(i + 1, sn):
                if line[i] * line[j] != 0:
                    shared[i][j] += 1
    line_total = np.sum(line_list, axis=0)
    return (line_total, shared)

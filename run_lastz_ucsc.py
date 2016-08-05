#!/usr/bin/env python
import os
import random
import string
import subprocess
import sys
#from itertools import zip
import ntpath
import sys
from time import sleep

def touch(fname):
    if os.path.exists(fname):
        os.utime(fname, None)
    else:
        open(fname, 'a').close()

def path_leaf(path):
    head, tail = ntpath.split(path)
    return head, tail

def exc_handler(tries_remaining, exception, delay):
    sys.stderr.write("Caught {}, {} tries remaining, sleeping for {} seconds".format(exception, tries_remaining, delay))

def retries(max_tries, delay=1, backoff=2, exceptions=(Exception,), hook=None):
    def dec(func):
        def f2(*args, **kwargs):
            mydelay = delay
            tries = range(max_tries)
            tries = list(reversed(tries))
            for tries_remaining in tries:
                try:
                   return func(*args, **kwargs)
                except exceptions as e:
                    if tries_remaining > 0:
                        if hook is not None:
                            hook(tries_remaining, e, mydelay)
                        sleep(mydelay)
                        mydelay = mydelay * backoff
                    else:
                        raise
                else:
                    break
        return f2
    return dec


TEMPLATE = """#PBS -S /bin/bash
#PBS -q cmb
#PBS -e {error_log}
#PBS -o {output_log}
#PBS -l vmem=10G
#PBS -l pmem=10G
#PBS -l mem=10G
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:30:00
#PBS -N {name}
set -eou pipefail
export PATH=/home/cmb-panasas2/skchoudh/software_frozen/anaconda2/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin
cd {WORK_DIR}
mkdir -p psl {tmpDir}
{BIN_DIR}/blastz-run-ucsc {TARGET} {QUERY} {DEF_FILE} psl/{out_filename}.psl.gz -outFormat psl -gz
"""

@retries(10, delay=60, hook=exc_handler)
def make_submission(job_command):
    output = subprocess.getoutput(job_command)
    if 'submit error' in output:
        raise Exception('Error with job {}'.format(job_command))
    return output
WORK_DIR='/home/cmb-panasas2/skchoudh/galGal4_vs_geoFor1_2bit_parallel_ucsc'
DEF_FILE='/home/cmb-panasas2/skchoudh/galGal4_vs_geoFor1_2bit_parallel_ucsc/DEF_FILE'
BIN_DIR='/home/cmb-panasas2/skchoudh/multizscripts/scripts'
TMP_DIR='/staging/as/skchoudh/scratch'
lastzParams='O=400 E=30 H=2000 L=2200 K=3000'




def walk_files(target_file, query_file):
    index1 = 0
    index2 = 0
    with open(target_file) as t:
        for index1, tline in enumerate(t):
            tline = tline.strip()
            with open(query_file) as q:
                for index2, qline in enumerate(q):
                    qline = qline.strip()
                    target = tline
                    query = qline
                    _, tname  = path_leaf(target)
                    _, qname = path_leaf(query)
                    job_name = '{}-{}'.format(index1, index2)

                    job_script = os.path.join(WORK_DIR, 'jobs', '{}-{}.sh'.format(index1, index2))
                    error_log = os.path.join(WORK_DIR, 'logs', '{}-{}.e'.format(index1, index2))
                    output_log = os.path.join(WORK_DIR, 'logs', '{}-{}.o'.format(index1, index2))

                    touch(error_log)
                    touch(output_log)
                    open(job_script, 'w').write(TEMPLATE.format(tmpDir=TMP_DIR,
                                                                WORK_DIR=WORK_DIR,
                                                                TARGET=target,
                                                                QUERY=query,
                                                                TNAME=tname,
                                                                QNAME=qname,
                                                                lastzParams=lastzParams,
                                                                error_log=error_log,
                                                                output_log=output_log,
                                                                out_filename='{}.{}'.format(tname, qname),
                                                                BIN_DIR=BIN_DIR,
                                                                DEF_FILE=DEF_FILE,
                                                                name='lastz_{}'.format(job_name)
                                                                ))
                    output = make_submission('qsub {}'.format(job_script))
                    print(output)

if __name__ == '__main__':
    target = sys.argv[1]
    query = sys.argv[2]
    walk_files(target, query)

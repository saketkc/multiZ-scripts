include:
    "config.py"

import errno
import glob
import os
import ntpath
import random
import string
import subprocess
from time import sleep

WORK_DIR = srcdir('')
BIN_DIR = os.path.abspath(BIN_DIR)

PREFIX_OUT = 'processed_data/{}_VS_{}'.format(TARGET, QUERY)
PREFIX_OUT_QUERY = PREFIX_OUT + '/' + QUERY
PREFIX_OUT_TARGET = PREFIX_OUT + '/' + TARGET

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def create_def_file(blastz_o=400,
                    blastz_e=30,
                    blastz_h=2000,
                    blastz_l=2200,
                    blastz_k=3000,
                    seq1_dir=None,
                    seq1_len=None,
                    seq2_dir=None,
                    seq2_len=None
                    ):
    def_file_template = """
    BLASTZ=lastz
    BLASTZ_O={blastz_o}
    BLASTZ_E={blastz_e}
    BLASTZ_H={blastz_h}
    BLASTZ_L={blastz_l}
    BLASTZ_K={blastz_k}

    SEQ1_DIR={seq1_dir}
    SEQ1_LEN={seq1_len}
    SEQ2_DIR={seq2_dir}
    SEQ2_LEN={seq2_len}

    TMPDIR=/staging/as/skchoudh
    """.format(blastz_o=blastz_o,
               blastz_e=blastz_e,
               blastz_h=blastz_h,
               blastz_l=blastz_l,
               blastz_k=blastz_k,
               seq1_dir=seq1_dir,
               seq1_len=seq1_len,
               seq2_dir=seq2_dir,
               seq2_len=seq2_len,
              )


def get_partlst_files(genome):
    return glob.glob('processed_data/{}_VS_{}/{}PartList/*.lst'.format(TARGET, QUERY, genome))

def get_partlst_filenames(genome):
    filelist = list(map(os.path.basename, get_partlst_files(genome)))
    filelist = [os.path.splitext(f)[0] for f in filelist]
    return filelist

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

@retries(20, delay=100, hook=exc_handler)
def make_submission(job_command):
    output = subprocess.getoutput(job_command)
    if 'submit error' in output:
        raise Exception('Error with job {}'.format(job_command))
    return output

rule all:
    input:
        expand('raw_data/{genome}.2bit', genome=GENOMES),
        expand('processed_data/{target}_VS_{query}/{genome}.chrom.sizes', genome=GENOMES, target=TARGET, query=QUERY),
        expand('processed_data/{target}_VS_{query}/{genome}.lst', genome=GENOMES, target=TARGET, query=QUERY),
        expand('processed_data/{target}_VS_{query}/{genome}PartList/', genome=TARGET, target=TARGET, query=QUERY),
        expand('processed_data/{target}_VS_{query}/{genome}PartList/', genome=QUERY, target=TARGET, query=QUERY),
        expand('processed_data/{target}_VS_{query}/{genome}PartList_2bit/', genome=TARGET, target=TARGET, query=QUERY),
        expand('processed_data/{target}_VS_{query}/{genome}PartList_2bit/', genome=QUERY, target=TARGET, query=QUERY),
        expand('processed_data/{target}_VS_{query}/psl', target=TARGET, query=QUERY),


rule install_requirements:
    shell: 'source activate {PYENV} & conda install -y {REQUIREMENTS}'

rule download_2bit:
    output: expand('raw_data/{genome}.2bit', genome=GENOMES)
    run:
        for genome in GENOMES:
            shell('wget -cP raw_data http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/{genome}.2bit')

rule create_chrominfo:
    input:
        expand('raw_data/{genome}.2bit', genome=GENOMES)
    output:
        expand('processed_data/{target}_VS_{query}/{genome}.chrom.sizes', genome=GENOMES, target=TARGET, query=QUERY),
    run:
        for index, genome in enumerate(GENOMES):
            input_f = input[index]
            output_f = output[index]
            shell('twoBitInfo {input_f} stdout | sort -k2nr > {output_f}')

rule create_target_partitions:
    input:
        'raw_data/{TARGET}.2bit',
        'processed_data/{TARGET}_VS_{QUERY}/{TARGET}.chrom.sizes'
    output:
        'processed_data/{TARGET}_VS_{QUERY}/{TARGET}.lst',
        'processed_data/{TARGET}_VS_{QUERY}/{TARGET}PartList/',
    shell:
        '{BIN_DIR}/partitionSequence.pl {TARGET_CHUNK} {TARGET_LAP} '
        '{input[0]} {input[1]} {TARGET_LIMIT} '
        '-lstDir {output[1]} > {output[0]}'

rule create_query_partitions:
    input:
        'raw_data/{QUERY}.2bit',
        'processed_data/{TARGET}_VS_{QUERY}/{QUERY}.chrom.sizes'
    output:
        'processed_data/{TARGET}_VS_{QUERY}/{QUERY}.lst',
        'processed_data/{TARGET}_VS_{QUERY}/{QUERY}PartList/',
    shell:
        '{BIN_DIR}/partitionSequence.pl {QUERY_CHUNK} {QUERY_LAP} '
        '{input[0]} {input[1]} {QUERY_LIMIT} '
        '-lstDir {output[1]} > {output[0]}'

rule create_target_lst_files:
    input:
        get_partlst_files(TARGET),
        'processed_data/{TARGET}_VS_{QUERY}/{TARGET}.lst'
    output:
        'processed_data/{TARGET}_VS_{QUERY}/{TARGET}PartList_2bit/',
    run:
        for tPart in input[:-1]:
            tPart = tPart.replace('.lst', '')
            shell('''sed -e 's#.*.2bit:##;' {tPart}.lst | twoBitToFa -seqList=stdin raw_data/{TARGET}.2bit stdout | faToTwoBit stdin {tPart}.2bit''')
            shell('''mkdir -p {output} && mv processed_data/{TARGET}_VS_{QUERY}/{TARGET}PartList/*.2bit {output}''')

rule create_query_lst_files:
    input:
        get_partlst_files(QUERY),
        'processed_data/{TARGET}_VS_{QUERY}/{QUERY}.lst'
    output:
        'processed_data/{TARGET}_VS_{QUERY}/{QUERY}PartList_2bit/',
    run:
        for qPart in input[:-1]:
            qPart = qPart.replace('.lst', '')
            shell('''sed -e 's#.*.2bit:##;' {qPart}.lst | twoBitToFa -seqList=stdin raw_data/{QUERY}.2bit stdout | faToTwoBit stdin {qPart}.2bit''')
            shell('''mkdir -p {output} && mv processed_data/{TARGET}_VS_{QUERY}/{QUERY}PartList/*.2bit {output}''')


rule create_psl:
    input:
        'processed_data/{TARGET}_VS_{QUERY}/{TARGET}.lst',
        'processed_data/{TARGET}_VS_{QUERY}/{QUERY}.lst',
        'processed_data/{TARGET}_VS_{QUERY}/{TARGET}PartList_2bit/',
        'processed_data/{TARGET}_VS_{QUERY}/{QUERY}PartList_2bit/',
    output:
        'processed_data/{TARGET}_VS_{QUERY}/psl',
    run:
        mkdir_p(os.path.join(WORK_DIR, 'jobs'))
        mkdir_p(os.path.join(WORK_DIR, 'logs'))

        with open(input[0]) as t:
            for index1, tline in enumerate(t):
                tline = tline.strip()
                with open(input[1]) as q:
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
                        open(job_script, 'w').write(BLASTZ_TEMPLATE.format(
                                                                    error_log=error_log,
                                                                    output_log=output_log,
                                                                    name='lastz_{}'.format(job_name),
                                                                    PATH=PATH,
                                                                    PYENV=PYENV,
                                                                    WORK_DIR=WORK_DIR,
                                                                    BIN_DIR=BIN_DIR,
                                                                    TARGET=target,
                                                                    QUERY=query,
                                                                    DEF_FILE=WORK_DIR+'/'+DEF_FILE,
                                                                    out_filename='{}.{}'.format(tname, qname),
                                                                    ))
                        output = make_submission('qsub {}'.format(job_script))

rule chainer:
    input:
        'processed_data/{TARGET}_VS_{QUERY}/{TARGET}.lst'
        'processed_data/{TARGET}_VS_{QUERY}/psl',
    output:
        'processed_data/{TARGET}_VS_{QUERY}/chain',
    run: #'''for T in `cat ${input[0]} | sed -e "s#${WORK_DIR}/##" | sed -e "s#${TARGET}PartList/##"`
        target_lines = open(input).readlines()
        for line in target_lines:
            line = line.strip()
            line = line.replace(WORK_DIR, '').replace('{}PartList/'.format(TARGET), '')
            job_name = '{}-{}'.format(line)
            job_script = os.path.join(WORK_DIR, 'jobs', '{}-{}.sh'.format(index1, index2))
            open(job_script, 'w').write(psl_fileprefix=line, 
                                        target=params['target_2bit'],
                                        query=params['query_2bit'])

            output = make_submission('qsub {}'.format(job_script))

rule clean_all:
    shell:
        'rm -rf processed_data jobs logs psl chain'

rule clean_data:
    shell:
        'rm -rf jobs logs psl chain'


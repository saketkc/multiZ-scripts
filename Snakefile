include:
    "config.py"

import glob
import os

WORK_DIR = srcdir('')
PREFIX_OUT = 'processed_data/{}_VS_{}'.format(TARGET, QUERY)
PREFIX_OUT_QUERY = PREFIX_OUT + '/' + QUERY
PREFIX_OUT_TARGET = PREFIX_OUT + '/' + TARGET


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

@retries(10, delay=60, hook=exc_handler)
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
    params:
        out_lst_file=expand('processed_data/{TARGET}_VS_{QUERY}/{TARGET}.lst', TARGET=TARGET, QUERY=QUERY),
        out_lst_dir=expand('processed_data/{TARGET}_VS_{QUERY}/{TARGET}PartList/', TARGET=TARGET, QUERY=QUERY),
    output:
        #'processed_data/{TARGET}_VS_{QUERY}/{TARGET}PartList/{partn}.lst',
        ##'processed_data/{TARGET}_VS_{QUERY}/{TARGET}PartList/',
        'processed_data/{TARGET}_VS_{QUERY}/{TARGET}.lst',
        'processed_data/{TARGET}_VS_{QUERY}/{TARGET}PartList/',
        #dynamic(PREFIX_OUT+'/'+TARGET+'PartList/{partn}.lst')
    shell:
        '{BIN_DIR}/partitionSequence.pl {TARGET_CHUNK} {TARGET_LAP} '
        '{input[0]} {input[1]} {TARGET_LIMIT} '
        '-lstDir {output[1]} > {output[0]}'

rule create_query_partitions:
    input:
        'raw_data/{QUERY}.2bit',
        'processed_data/{TARGET}_VS_{QUERY}/{QUERY}.chrom.sizes'
    params:
        out_lst_file=expand('processed_data/{TARGET}_VS_{QUERY}/{QUERY}.lst', TARGET=TARGET, QUERY=QUERY),
        out_lst_dir=expand('processed_data/{TARGET}_VS_{QUERY}/{QUERY}PartList/', TARGET=TARGET, QUERY=QUERY)
    output:
        #'processed_data/{TARGET}_VS_{QUERY}/{QUERY}PartList/{partn}.lst',
        'processed_data/{TARGET}_VS_{QUERY}/{QUERY}.lst',
        'processed_data/{TARGET}_VS_{QUERY}/{QUERY}PartList/',
        #PREFIX_OUT_QUERY+'PartList/{part,}.lst'
        #dynamic(PREFIX_OUT+'/'+QUERY+'PartList/{partn}.lst')
    shell:
        '{BIN_DIR}/partitionSequence.pl {QUERY_CHUNK} {QUERY_LAP} '
        '{input[0]} {input[1]} {QUERY_LIMIT} '
        '-lstDir {output[1]} > {output[0]}'

rule create_target_lst_files:
    input:
        get_partlst_files(TARGET),
        'processed_data/{TARGET}_VS_{QUERY}/{TARGET}.lst'
        #dynamic(PREFIX_OUT+'/'+TARGET+'PartList/{partn}.lst')
    params:
        output=expand('processed_data/{TARGET}_VS_{QUERY}/{TARGET}.partlst', TARGET=TARGET, QUERY=QUERY)
    output:
        #dynamic(PREFIX_OUT+'/'+TARGET+'PartList/{partn}.2bit')
        'processed_data/{TARGET}_VS_{QUERY}/{TARGET}PartList_2bit/',
        ##expand('processed_data/{TARGET}_VS_{QUERY}/{TARGET}.partlst', TARGET=TARGET, QUERY=QUERY),
        #dynamic(PREFIX_OUT+'/'+TARGET+'PartList/{partn}.2bit')
        ##expand('processed_data/{target}_VS_{query}/{query}PartList/{sample}.2bit',
        ##       query=QUERY,
        ##      target=TARGET,
        ##       sample=get_partlst_filenames(TARGET))

    run:
        for tPart in input[:-1]:
            #print(tPart)
            tPart = tPart.replace('.lst', '')
            shell('''sed -e 's#.*.2bit:##;' {tPart}.lst | twoBitToFa -seqList=stdin raw_data/{TARGET}.2bit stdout | faToTwoBit stdin {tPart}.2bit''')
            shell('''mkdir -p processed_data/{TARGET}_VS_{QUERY}/{TARGET}PartList_2bit && mv processed_data/{TARGET}_VS_{QUERY}/{TARGET}PartList/*.2bit {TARGET}PartList_2bit''')

rule create_query_lst_files:
    priority: 1
    input:
        get_partlst_files(QUERY),
        'processed_data/{TARGET}_VS_{QUERY}/{QUERY}.lst'
    output:
        'processed_data/{TARGET}_VS_{QUERY}/{QUERY}PartList_2bit/',
    run:
        for qPart in input[:-1]:
            qPart = qPart.replace('.lst', '')
            shell('''sed -e 's#.*.2bit:##;' {qPart}.lst | twoBitToFa -seqList=stdin raw_data/{QUERY}.2bit stdout | faToTwoBit stdin {qPart}.2bit''')
            shell('''mkdir -p processed_data/{TARGET}_VS_{QUERY}/{QUERY}PartList_2bit && mv processed_data/{TARGET}_VS_{QUERY}/{QUERY}PartList/*.2bit {TARGET}PartList_2bit''')
rule clean:
    shell:
        'rm -rf processed_data'

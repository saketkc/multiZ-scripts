include:
    "config.py"

import glob
import os

WORK_DIR = srcdir('')

def get_partlst_files(genome):
    return glob.glob('processed_data/{}/{}PartList/*.lst'.format(genome, genome))

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
        expand('processed_data/{target}VS{query}/{genome}.chrom.sizes', genome=GENOMES, target=TARGET, query=QUERY),
        expand('processed_data/{target}VS{query}/{genome}.lst', genome=GENOMES, target=TARGET, query=QUERY),
        expand('processed_data/{target}VS{query}/{genome}PartList/part000.2bit',  genome=QUERY, target=TARGET, query=QUERY),
        expand('processed_data/{target}VS{query}/{genome}PartList/part000.2bit',  genome=TARGET, target=TARGET, query=QUERY),

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
        expand('processed_data/{target}VS{query}/{genome}.chrom.sizes', genome=GENOMES, target=TARGET, query=QUERY),
    run:
        for index, genome in enumerate(GENOMES):
            input_f = input[index]
            output_f = output[index]
            shell('twoBitInfo {input_f} stdout | sort -k2nr > {output_f}')

rule create_target_partitions:
    input:
        'raw_data/{TARGET}.2bit',
        'processed_data/{TARGET}VS{QUERY}/{TARGET}.chrom.sizes'
    output:
        'processed_data/{TARGET}VS{QUERY}/{TARGET}PartList/',
        'processed_data/{TARGET}VS{QUERY}/{TARGET}.lst'
    shell:
        '{BINDIR}/partitionSequence.pl {TARGET_CHUNK} {TARGET_LAP} '
        '{input[0]} {input[1]} {TARGET_LIMIT} '
        '-lstDir {output[0]} > {output[1]}'

rule create_query_partitions:
    input:
        'raw_data/{QUERY}.2bit',
        'processed_data/{TARGET}VS{QUERY}/{QUERY}.chrom.sizes'
    output:
        'processed_data/{TARGET}VS{QUERY}/{QUERY}PartList/',
        'processed_data/{TARGET}VS{QUERY}/{QUERY}.lst'
    shell:
        '{BINDIR}/partitionSequence.pl {QUERY_CHUNK} {QUERY_LAP} '
        '{input[0]} {input[1]} {QUERY_LIMIT} '
        '-lstDir {output[0]} > {output[1]}'

rule create_query_lst_files:
    input: get_partlst_files(QUERY)
    params:
        query='processed_data/{TARGET}VS{QUERY}/{QUERY}PartList',
    output:
        expand('processed_data/{target}VS{query}/{query}PartList/{sample}.2bit',
               query=QUERY,
               target=TARGET,
               sample=get_partlst_filenames(QUERY))
    run:
        for qPart in input:
            qPart = qPart.replace('.lst', '')
            shell('''sed -e 's#.*.2bit:##;' {qPart}.lst | twoBitToFa -seqList=stdin raw_data/{QUERY}.2bit stdout | faToTwoBit stdin {qPart}.2bit''')

rule create_target_lst_files:
    input: get_partlst_files(TARGET)
    params:
        target='processed_data/{TARGET}VS{QUERY}/{TARGET}PartList',
    output:
        expand('processed_data/{target}VS{query}/{target}PartList/{sample}.2bit',
               query=QUERY,
               target=TARGET,
               sample=get_partlst_filenames(TARGET))
    run:
        for tPart in input:
            tPart = tPart.replace('.lst', '')
            shell('''sed -e 's#.*.2bit:##;' {tPart}.lst | twoBitToFa -seqList=stdin raw_data/{TARGET}.2bit stdout | faToTwoBit stdin {tPart}.2bit''')

rule create_psl:
    input:
        'processed_data/{TARGET}VS{QUERY}/{TARGET}.lst',
        'processed_data/{TARGET}VS{QUERY}/{QUERY}.lst'
    output:
        'processed_data/{TARGET}VS{QUERY}/psl',
    run:
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
                        open(job_script, 'w').write(BLASTZ_TEMPLATE.format(tmpDir=TMP_DIR,
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

rule chainer:
    input:
        'processed_data/{TARGET}VS{QUERY}/{TARGET}.lst'
    params: 
        target_2bit='raw_data/{TARGET}.2bit',
        query_2bit='raw_data/{QUERY}.2bit',
    output:
        dynamic('processed_data/'+TARGET'vs'+QUERY+'/chain/')
    run: #'''for T in `cat ${input[0]} | sed -e "s#${WORKDIR}/##" | sed -e "s#galGal4PartList/##"`
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
rule clean:
    shell:
        'rm -rf processed_data'

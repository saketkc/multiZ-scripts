import glob, os

BINDIR = 'bin'
TARGET = 'galGal4'
QUERY = 'mm10'

GENOMES = [TARGET, QUERY]
TARGET_CHUNK = '10000000'
TARGET_LAP = '10000'
TARGET_LIMIT = '300'

QUERY_CHUNK = '20000000'
QUERY_LAP = '0'
QUERY_LIMIT = '300'

rule all:
    input: 
        expand('raw_data/{genome}.2bit', genome=GENOMES),
        expand('processed_data/{genome}/{genome}.chrom.sizes', genome=GENOMES),
        expand('processed_data/{genome}/{genome}.lst', genome=GENOMES),
        expand('processed_data/{genome}/{genome}PartList/part000.2bit',  genome=QUERY),
        expand('processed_data/{genome}/{genome}PartList/part000.2bit',  genome=TARGET)
        #expand('processed_data/{genome}PartList', genome=GENOMES)

rule download_2bit:
    output: expand('raw_data/{genome}.2bit', genome=GENOMES)
    run:
        for genome in GENOMES:
            shell('wget -cP raw_data http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/bigZips/{genome}.2bit')

rule create_chrominfo:
    input: expand('raw_data/{genome}.2bit', genome=GENOMES)
    output: expand('processed_data/{genome}/{genome}.chrom.sizes', genome=GENOMES)
    run:
        for index, genome in enumerate(GENOMES):
            input_f = input[index]
            output_f = output[index]
            shell('twoBitInfo {input_f} stdout | sort -k2nr > {output_f}')

rule create_target_partitions:
    input:
        expand('raw_data/{target}.2bit', target=TARGET),
        expand('processed_data/{target}/{target}.chrom.sizes', target=TARGET),
    output:
        expand('processed_data/{target}/{target}PartList', target=TARGET),
        expand('processed_data/{target}/{target}.lst', target=TARGET)
    shell:
        '{BINDIR}/partitionSequence.pl {TARGET_CHUNK} {TARGET_LAP} '
        '{input[0]} {input[1]} {TARGET_LIMIT} '
        '-lstDir {output[0]} > {output[1]}'

rule create_query_partitions:
    input:
        expand('raw_data/{query}.2bit', query=QUERY),
        expand('processed_data/{query}/{query}.chrom.sizes', query=QUERY),
    output:
        expand('processed_data/{query}/{query}PartList', query=QUERY),
        expand('processed_data/{query}/{query}.lst', query=QUERY)
    shell:
        '{BINDIR}/partitionSequence.pl {QUERY_CHUNK} {QUERY_LAP} '
        '{input[0]} {input[1]} {QUERY_LIMIT} '
        '-lstDir {output[0]} > {output[1]}'


def get_partlst_files(genome):
    return glob.glob('processed_data/{}/{}PartList/*.lst'.format(genome, genome))

def get_partlst_filenames(genome):
    filelist = list(map(os.path.basename, get_partlst_files(genome)))
    filelist = [os.path.splitext(f)[0] for f in filelist]
    return filelist

rule create_query_lst_files:
    input: get_partlst_files(QUERY)
    params:
        query='processed_data/{}/{}PartList'.format(QUERY, QUERY),
    output:
        expand('processed_data/{query}/{query}PartList/{sample}.2bit',
               query=QUERY,
               sample=get_partlst_filenames(QUERY))
    run:
        for qPart in input:
            qPart = qPart.replace('.lst', '')
            shell('''sed -e 's#.*.2bit:##;' {qPart}.lst | twoBitToFa -seqList=stdin raw_data/{QUERY}.2bit stdout | faToTwoBit stdin {qPart}.2bit''')

rule create_target_lst_files:
    input: get_partlst_files(TARGET)
    params:
        target='processed_data/{}/{}PartList'.format(TARGET, TARGET)
    output:
        expand('processed_data/{target}/{target}PartList/{sample}.2bit',
               target=TARGET,
               sample=get_partlst_filenames(TARGET))
    run:
        for tPart in input:
            tPart = tPart.replace('.lst', '')
            shell('''sed -e 's#.*.2bit:##;' {tPart}.lst | twoBitToFa -seqList=stdin raw_data/{TARGET}.2bit stdout | faToTwoBit stdin {tPart}.2bit''')

rule clean:
    shell:
        'rm -rf processed_data'

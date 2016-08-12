PYENV = 'py35'

PATH = '/home/cmb-panasas2/skchoudh/software_frozen/anaconda2/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin'

DEF_FILE = './DEF_FILE'

TARGET = 'galGal4'
QUERY = 'geoFor1'

BIN_DIR = './bin'

GENOMES = [TARGET, QUERY]

TARGET_CHUNK = '10000000'
TARGET_LAP = '10000'
TARGET_LIMIT = '300'

QUERY_CHUNK = '20000000'
QUERY_LAP = '0'
QUERY_LIMIT = '300'

REQUIREMENTS = ['ucsc-axtchain',
                'ucsc-axttomaf',
                'ucsc-axttopsl',
                'ucsc-chainantirepeat',
                'ucsc-chainmergesort',
                'ucsc-chainnet',
                'ucsc-chainprenet',
                'ucsc-fasplit',
                'ucsc-lavtopsl',
                'ucsc-liftup',
                'ucsc-netchainsubset',
                'ucsc-netclass',
                'ucsc-nettoaxt',
                'ucsc-pslcheck',
                'ucsc-twobitinfo',
                'ucsc-twobittofa']

TEMPLATE_HEADER = """#PBS -S /bin/bash
#PBS -q cmb
#PBS -e {error_log}
#PBS -o {output_log}
#PBS -l vmem=10G
#PBS -l pmem=10G
#PBS -l mem=10G
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:30:00
#PBS -N {name}
#set -eou pipefail
export PATH={PATH}
source activate {PYENV}
cd {WORK_DIR}"""

BLASTZ_TEMPLATE = TEMPLATE_HEADER + """
mkdir -p psl
{BIN_DIR}/blastz-run-ucsc {TARGET} {QUERY} {DEF_FILE} psl/{out_filename}.psl.gz -outFormat psl -gz
"""
#.format(BIN_DIR=BIN_DIR, PYENV=PYENV)

CHAIN_TEMPLATE = TEMPLATE_HEADER +"""
zcat psl/{psl_fileprefix}.*.psl.gz \
| axtChain -psl -verbose=0 {chainParams} stdin {TARGET} {QUERY} stdout \
| chainAntiRepeat {TARGET} {QUERY} stdin chain/{psl_fileprefix}.chain
"""

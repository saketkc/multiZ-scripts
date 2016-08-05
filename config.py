PYENV = 'py35'

TARGET = 'galGal4'
QUERY = 'geoFor1'

BINDIR = './bin'

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

import os
from genomes import query_genomes, target_genome
import errno
import shutil
import re

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            print ('....')
            raise

for query in query_genomes:

    dest_dir = os.path.join('../', query)
    source_dir = os.path.abspath('../config_template')
    config_temp = os.path.join(source_dir, 'config.py')
    mkfile_temp = os.path.join('../', 'Snakefile')

    with open(config_temp, 'r') as f:
        config_text = str(f.read().strip())
    config_text = config_text.replace('QUERY =', 'QUERY = \'{}\''.format(query))
    config_text = config_text.replace('TARGET =', 'TARGET = \'{}\''.format(target_genome))

    mkdir_p(dest_dir)
    shutil.copy(mkfile_temp, dest_dir)
    mkdir_p(os.path.join(dest_dir, 'config'))
    with open(os.path.join(dest_dir, 'config', 'config.py'), 'w') as f:
        f.write(config_text)



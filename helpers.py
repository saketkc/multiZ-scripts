import glob
import os

def get_partlst_files(genome):
    return glob.glob('processed_data/{}/{}PartList/*.lst'.format(genome, genome))

def get_partlst_filenames(genome):
    filelist = list(map(os.path.basename, get_partlst_files(genome)))
    filelist = [os.path.splitext(f)[0] for f in filelist]
    return filelist

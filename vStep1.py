import multiprocessing
from functools import partial
from vFunctions import set_variables
from vFunctions import align_reads
import glob
import os

os.chdir('/Users/tstuber/Desktop/fastq_test/test')

def read_aligner(species_call, single_directory):
    os.chdir(single_directory)
    R1 = glob.glob('*_R1*fastq.gz')[0]
    R2 = glob.glob('*_R2*fastq.gz')[0]
    print("R1 and R2: %s %s" % (R1, R2))
    try:
        if species_call:
            sample_attributes = set_variables(R1, R2, species_call) #4
        else:
            species_call = False
            sample_attributes = set_variables(R1, R2, species_call) #4
    except:
        species_call = False
        sample_attributes = set_variables(R1, R2, species_call) #4
    stat_summary = align_reads(sample_attributes) #5
    return stat_summary

def pass_var():
    run_list = glob.glob('*') #1

    # with futures.ProcessPoolExecutor(max_workers=limited_cpu_count) as pool: 
    #     for stat_summary in pool.map(read_aligner, run_list):
    #         master_stat_summary.append(stat_summary)

    master_stat_summary = []
    species_call = 'af'
    pool = multiprocessing.Pool()
    func = partial(read_aligner, species_call)
    for stat_summary in pool.map(func, run_list):
        master_stat_summary.append(stat_summary)

    pool.close()
    pool.join()

pass_var()
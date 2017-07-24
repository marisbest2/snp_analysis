#!/usr/bin/env python

import glob

R1 = glob.glob('*_R1*fastq.gz')
R2 = glob.glob('*_R2*fastq.gz')
list_len = len(list_of_files)
if (list_len % 2 != 0):
    print("\n#####Check paired files.  Unpaired files seen by odd number of counted FASTQs\n\n")
    sys.exit(0)



fastqs = glob.glob(zips + '/*.vcf')

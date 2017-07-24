#!/usr/bin/env python

import glob
import sys

all_file_types_count = len(glob.glob('*.*'))
fastq_check = len(glob.glob('*fastq.gz'))
vcf_check = len(glob.glob('*vcf'))

if fastq_check > 0:
    fastq_check = True
if vcf_check > 0:
    vcf_check = True
if fastq_check and vcf_check:
    print("\n#####You have a mix of FASTQ and VCF files.  This is not allowed\n\n")
    sys.exit(0)

if fastq_check:
    R1 = glob.glob('*_R1*fastq.gz')
    R2 = glob.glob('*_R2*fastq.gz')
    R1count = len(R1)
    R2count = len(R2)
    fastq_count = R1count + R2count
    if (fastq_count % 2 != 0):
        print("\n#####Check paired files.  Unpaired files seen by odd number of counted FASTQs\n\n")
        sys.exit(0)
    if (R1count != R2count):
        print("\n#####Check paired files.  R1 files do not equal R2\n\n")
        sys.exit(0)
    if (all_file_types_count != fastq_count):
        print("\n#####Only zipped FASTQ files are allowed in directory\n\n")
        sys.exit(0)
    elif (fastq_count > 1):
        print("execute script 1")
elif vcf_check:
    vcfs_count = len(glob.glob('*vcf'))
    if (all_file_types_count != vcfs_count):
        print("\n#####You have more than just VCF files in your directory.  Only VCF files are allowed if running script 2\n\n")
        sys.exit(0)
else:
    print ("######There was an error determining the file type.")
    sys.exit(0)



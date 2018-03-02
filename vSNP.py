#!/usr/bin/env python

import os
import sys
import argparse
import glob
import json
import multiprocessing
from concurrent import futures
import textwrap
from datetime import datetime

from vFunctions import make_global
from vFunctions import get_species
from vFunctions import send_email
from vFunctions import run_script2
from vFunctions import loop_all
from vFunctions import loop_resticted
from vFunctions import read_aligner
from vFunctions import make_species_call_global

startTime = datetime.now()
print ("\n\n*** START ***\n")
print ("Start time: %s" % startTime)

home = os.path.expanduser("~")
root_dir = str(os.getcwd())
cpu_count = multiprocessing.cpu_count()
limited_cpu_count = int(cpu_count/6)
if limited_cpu_count == 0:
    limited_cpu_count = 1

parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

---------------------------------------------------------
vSNP --> get SNPs, group SNPs, verify SNPs

vSNP is called on a working directory containing FASTQ or VCF files.

See documentation at: https://usda-vs.github.io/snp_analysis/

        Step 1: FASTQs --> VCF

        Step 2: VCFs --> Tables & Trees

-s <OPTIONAL SPECIES TYPES>: af, h37, ab1, ab3, suis1, mel1, mel1b, mel2, mel3, canis, ceti1, ceti2, ovis, neo, para, salmonella

'''), epilog='''---------------------------------------------------------''')

#universal
parser.add_argument('-s', '--species', action='store', dest='species', help='OPTIONAL: Used to FORCE species type <see options above>')

parser.add_argument('-d', '--debug', action='store_true', dest='debug_call', help='debug, run without loop.map for loops')
parser.add_argument('-a', '--all_vcf', action='store_true', dest='all_vcf', help='make tree using all VCFs')
parser.add_argument('-e', '--elite', action='store_true', dest='elite', help='create a tree with on elite sample representation')
parser.add_argument('-f', '--filter', action='store_true', dest='filter', help='Find possible positions to filter')

parser.add_argument('-q', '--quiet', action='store_true', dest='quiet', help='[**APHIS only**] prevent stats going to cumlative collection')
parser.add_argument('-m', '--email', action='store', dest='email', help='[**APHIS only**, specify own SMTP address for functionality] email options: all, s, tod, jess, suelee, chris, email_address')
parser.add_argument('-u', '--upload', action='store_true', dest='upload', help='[**APHIS only**, specify own storage for functionality] upload files to the bioinfo drive')

parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.0.1')

args = parser.parse_args()
species_call = args.species
debug_call = args.debug_call
quiet_call = args.quiet
print ("\nSET ARGUMENTS: ")
print (args)
print("")

if args.email == "all":
    email_list = "tod.p.stuber@aphis.usda.gov, Jessica.A.Hicks@aphis.usda.gov, Christine.R.Quance@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov, patrick.m.camp@aphis.usda.gov, David.T.Farrell@aphis.usda.gov, Robin.L.Swanson@aphis.usda.gov, hannah.m.tharp@aphis.usda.gov, Doris.M.Bravo@aphis.usda.gov, eto3@cdc.gov"
elif args.email == "tod":
    email_list = "tod.p.stuber@aphis.usda.gov"
elif args.email == "jess":
    email_list = "Jessica.A.Hicks@aphis.usda.gov"
elif args.email == "suelee":
    email_list = "tod.p.stuber@aphis.usda.gov, Jessica.A.Hicks@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov, Doris.M.Bravo@aphis.usda.gov"
elif args.email == "chris":
    email_list = "tod.p.stuber@aphis.usda.gov, Jessica.A.Hicks@aphis.usda.gov, Christine.R.Quance@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov, eto3@cdc.gov"
elif args.email == "doris":
    email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, doris.m.bravo@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"
else:
    email_list = "tod.p.stuber@aphis.usda.gov"

##############################

fastq_check = 0
vcf_check = 0
all_file_types_count = len(glob.glob('*.*'))
fastq_check = len(glob.glob('*fastq.gz'))
vcf_check = len(glob.glob('*vcf'))

make_global(home, startTime, root_dir, cpu_count, limited_cpu_count, debug_call, quiet_call)

if species_call:
    make_species_call_global(species_call)

master_stat_summary = []
if fastq_check > 0 and vcf_check == 0:
    list_of_files = glob.glob('*gz')
    directory_list = loop_all(list_of_files) #1
    run_list = loop_resticted(directory_list) #2
    if debug_call:
        for single_directory in run_list:
            print("DEBUGGING, SAMPLES RAN INDIVIDUALLY")
            stat_summary = read_aligner(single_directory) #3
            master_stat_summary.append(stat_summary)
    else:
        print("SAMPLES RAN IN PARALLEL")
        with futures.ProcessPoolExecutor(max_workers=limited_cpu_count) as pool: 
            for stat_summary in pool.map(read_aligner, run_list):
                master_stat_summary.append(stat_summary)
    print("Done")
    with open("master_stat_summary.json", 'w') as outfile:
        json.dump(master_stat_summary, outfile)
    outfile.close()

elif vcf_check > 0 and fastq_check == 0:
    if not species_call:
        print("Checking VCF file for species type...")
        species_call = get_species()
        print("species_call %s" % species_call)
    vcfs_count = len(glob.glob('*vcf'))
    if (all_file_types_count != vcfs_count):
        print("\n#####You have more than just VCF files in your directory.  Only VCF files are allowed if running script 2\n\n")
        sys.exit(0)
    else:
        if quiet_call:
            print("#####Incorrect use of options when running script 2")
            sys.exit(0)
        else:
            if species_call:
                print("\n--> RUNNING SCRIPT 2\n") #
                run_script2(root_dir)
            else:
                print("#####Based on VCF CHROM id (reference used to build VCF) a matching species cannot be found neither was there a -s option given")
                sys.exit(0)
else:
    print("\n#####You have a mix of FASTQ and VCF files.  This is not allowed\n\n")
    sys.exit(0)

if args.email:
    send_email(email_list)

runtime = (datetime.now() - startTime)
print ("\n\nruntime: %s:  \n" % runtime)
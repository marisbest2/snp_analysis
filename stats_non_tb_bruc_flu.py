#!/usr/bin/env python

import os
import sys
import subprocess
import glob
import re
import shutil
import time
from datetime import datetime
from collections import OrderedDict
import csv
import xlrd
import vcf
import pandas as pd # pandas 0.18.1 tested, 0.19.2 does not work
import numpy as np
import xlsxwriter

access=False
root = str(os.getcwd())

def sizeof_fmt(num, suffix='B'):
    for unit in ['','K','M','G','T','P','E','Z']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)

ts = time.time()
st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')

#Cumulative stats
path_found = False
if os.path.isdir("/bioinfo11/TStuber/Results/stats"): #check bioinfo from server
    path_found = True
    copy_to = "/bioinfo11/TStuber/Results/stats"
elif os.path.isdir("/Volumes/root/TStuber/Results"): #check bioinfo from Mac
    path_found = True
    copy_to = "/Volumes/root/TStuber/Results/stats"
if path_found:
    summary_cumulative_file = copy_to + '/stat_alignment_culmulative_summary' + '.xlsx'
else:
    print("Bioinfo not connected")
    sys.exit(0)

###
summary_cumulative_file_temp = copy_to + '/stat_alignment_culmulative_summary' + st + '.xlsx'
workbook = xlsxwriter.Workbook(summary_cumulative_file_temp)
worksheet = workbook.add_worksheet()

row = 0
col = 0

top_row_header = ["time_stamp", "sample_name", "self.species", "reference_sequence_name", "R1size", "R2size", "allbam_mapped_reads", "genome_coverage", "ave_coverage", "ave_read_length", "unmapped_reads", "unmapped_assembled_contigs", "good_snp_count", "mlst_type", "octalcode", "sbcode", "hexadecimal_code", "binarycode"]
for header in top_row_header:
    worksheet.write(row, col, header)
    col += 1
###

list_of_files = glob.glob('*gz')
list_len = len(list_of_files)
if (list_len % 2 != 0):
    print("\n#####Check paired files.  Unpaired files seen by odd number of counted FASTQs\n\n")
    sys.exit(0)

for file in list_of_files:
    prefix_name=re.sub('_.*', '', file)
    print(prefix_name)
    if not os.path.exists(prefix_name):
        os.makedirs(prefix_name)
    shutil.move(file, prefix_name)

directory_list=[]
for f in  os.listdir('.'):
    if not f.startswith('.'):
        directory_list.append(f)

for d in directory_list:
    R1 = glob.glob(root + "/" + d + '/*_R1*fastq.gz')[0] #should just be one found.  puts the single list element to string
    R2 = glob.glob(root + "/" + d + '/*_R2*fastq.gz')[0]

    read_base = os.path.basename(R1)
    sample_name=re.sub('_.*', '', read_base)
    R1size = sizeof_fmt(os.path.getsize(R1))
    R2size = sizeof_fmt(os.path.getsize(R2))

    stat_summary = {}
    stat_summary["01"] = st
    stat_summary["02"] = sample_name
    stat_summary["03"] = "NA"
    stat_summary["04"] = "NA"
    stat_summary["05"] = R1size
    stat_summary["06"] = R2size
    stat_summary["07"] = "NON TB - BRUC SAMPLE *****************************************"

    stat_summary=OrderedDict(sorted(stat_summary.items()))
    for k, v in stat_summary.items():
        print("k: %s, v: %s" % (k, v))

    mytable=pd.read_excel(summary_cumulative_file)
    top=mytable[0:0]
    stat_summary_list = list(stat_summary.values())
    stat_summary_list.insert(0, st)
    while len(stat_summary_list) < 18:
        stat_summary_list.append("-")
    in_dict_stats=OrderedDict(zip(top, stat_summary_list))
    df_dict_stats=pd.DataFrame(in_dict_stats, index=[0])
    frames = [mytable, df_dict_stats]
    df_cum_stat = pd.concat(frames)
    #write over previous
    try:
        df_cum_stat.to_excel(summary_cumulative_file)
        access=True
    except:
        col = 0
        row += 1
        #run stats
        for v in stat_summary.values():
            worksheet.write(row, col, v)
            col += 1
workbook.close()

# if there was access to the cumlative stat results this temp can be deleted
if access:
    os.remove(summary_cumulative_file_temp)

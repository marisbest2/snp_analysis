#!/usr/bin/env python

import os
import sys
import glob
import re
import shutil
import time
import pandas as pd
from collections import OrderedDict
from datetime import datetime
from optparse import OptionParser
from concurrent import futures
import multiprocessing
import xlsxwriter
import script1
import smtplib,ssl
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email.mime.text import MIMEText
from email.utils import formatdate
from email import encoders

home = os.path.expanduser("~")

root = str(os.getcwd())
cpu_count = multiprocessing.cpu_count()

startTime = datetime.now()
print ("\n\n*** START ***\n")
print ("Start time: %s" % startTime)

###############################################

parser = OptionParser()
parser.add_option('-s', '--species', action='store', dest='species', help='--> -s option: USE TO FORCE SPECIES TYPE <--', metavar='<OPTIONAL options: bovis, h37, ab1, ab3, suis1, mel1, mel2, mel3, canis, ceti1, ceti2, para')
parser.add_option('-d', '--debug', action='store_true', dest='debug_call', help='debug, run without loop')
parser.add_option('-q', '--quiet', action='store_true', dest='quiet', help='prevent stats going to cumlative collection')
parser.add_option('-m', '--email', action='store', dest='email', help='email recipients: all, tod, jess, suelee, chris')

(options, args) = parser.parse_args()
print ("SET OPTIONS: ")
print (options)

if options.email == "all":
    email_list = "tod.p.stuber@aphis.usda.gov, Jessica.A.Hicks@aphis.usda.gov, Christine.R.Quance@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov, patrick.m.camp@aphis.usda.gov, David.T.Farrell@aphis.usda.gov, Robin.L.Swanson@aphis.usda.gov, hannah.m.tharp@aphis.usda.gov, Doris.M.Bravo@aphis.usda.gov"
elif options.email == "tod":
    email_list = "tod.p.stuber@aphis.usda.gov"
elif options.email == "jess":
    email_list = "Jessica.A.Hicks@aphis.usda.gov"
elif options.email == "suelee":
    email_list = "tod.p.stuber@aphis.usda.gov, Jessica.A.Hicks@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov"
elif options.email == "chris":
    email_list = "tod.p.stuber@aphis.usda.gov, Jessica.A.Hicks@aphis.usda.gov, Christine.R.Quance@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov"
else:
    email_list = "tod.p.stuber@aphis.usda.gov"

###############################################

list_of_files = glob.glob('*gz')
list_len = len(list_of_files)
if (list_len % 2 != 0):
    print("\n#####Check paired files.  Unpaired files seen by odd number of counted FASTQs\n\n")
    sys.exit(0)

sample_restriction = int(cpu_count/3)

for file in list_of_files:
    prefix_name=re.sub('_.*', '', file)
    print(prefix_name)
    if not os.path.exists(prefix_name):
        os.makedirs(prefix_name)
    shutil.move(file, prefix_name)

###
#Run stats
ts = time.time()
st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')
summary_file = root + '/stat_alignment_summary_' + st + '.xlsx'

workbook = xlsxwriter.Workbook(summary_file)
worksheet = workbook.add_worksheet()

row = 0
col = 0

top_row_header = ["sample_name", "self.species", "reference_sequence_name", "R1size", "R2size", "allbam_mapped_reads", "genome_coverage", "ave_coverage", "ave_read_length", "unmapped_reads", "unmapped_assembled_contigs", "good_snp_count", "mlst_type", "octalcode", "sbcode", "hexadecimal_code", "binarycode"]

for header in top_row_header:
    worksheet.write(row, col, header)
    col += 1
###
#Cumulative stats
path_found = False
if os.path.isdir("/bioinfo11/TStuber/Results/stats"): #check bioinfo from server
    path_found = True
    copy_to = "/bioinfo11/TStuber/Results/stats"
elif os.path.isdir("/Volumes/root/TStuber/Results"): #check bioinfo from Mac
    path_found = True
    copy_to = "/Volumes/root/TStuber/Results/stats"
else:
    copy_to="no_path"
    print("Bioinfo not connected")

if path_found:
    summary_cumulative_file = copy_to + '/stat_alignment_culmulative_summary' + '.xlsx'
###

def read_aligner(directory):
    os.chdir(directory)
    R1 = glob.glob('*_R1*fastq.gz')
    R2 = glob.glob('*_R2*fastq.gz')
    if options.species:
        sample = script1.script1(R1[0], R2[0], options.species) #force species
    else:
        sample = script1.script1(R1[0], R2[0]) #no species give, will find best
    stat_summary = sample.align_reads()
    return(stat_summary)

directory_list=[]
for f in  os.listdir('.'):
    if not f.startswith('.'):
        directory_list.append(f)

total_samples = len(directory_list)
sample_restriction = int(cpu_count/4)
lower_count = 0
upper_count = 1
while lower_count < total_samples:
    upper_count = lower_count + sample_restriction
    run_list = directory_list[lower_count:upper_count] #create a run list
    for i in run_list:
        directory_list.remove(i)
    total_samples = len(directory_list)
    print(run_list)

    print("Iterating directories")
    if options.debug_call: #run just one sample at a time to debug
        for d in run_list:
            print("DEBUGGING, SAMPLES RAN INDIVIDUALLY")
            stat_summary = read_aligner(d)
            col = 0
            row += 1
            for v in stat_summary.values():
                worksheet.write(row, col, v)
                col += 1
            os.chdir(root)
    else: # run all in run_list in parallel
        print("SAMPLES RAN IN PARALLEL")
        with futures.ProcessPoolExecutor(max_workers=sample_restriction) as pool: #max_workers=cpu_count
            for stat_summary in pool.map(read_aligner, run_list): #run in parallel run_list in read_aligner (script1)
                print("statsummary")
                print(stat_summary) #stat_summary returned from script 1
                col = 0
                row += 1
                #run stats
                for v in stat_summary.values():
                    worksheet.write(row, col, v) #stat summary to be attached in email and left in working directory
                    col += 1
                stats_lock = False
                if not options.quiet and path_found:
                    try: #try to open cumulative stats file
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
                        df_cum_stat.to_excel(summary_cumulative_file)
                    except: # if cannot open cumulative stats file place a copy of run stats in the stats folder
                        print("Cumulative stats file is locked")
                        stats_lock = True

try:
    runtime = (datetime.now() - startTime)
    col = 0
    row += 1
    worksheet.write(row, col, "runtime: %s: " % runtime)
    workbook.close()
except:
    print("ERROR CLOSING STATS FILE")
    pass

if stats_lock: #if file was locked try to copy
    try: #will not copy if path unavailable
        shutil.copy2(summary_file, copy_to)
    except:
        pass

####send email:

def send_email(email_list):
    text = "See attached:  "
    send_from = "tod.p.stuber@aphis.usda.gov"
    send_to = email_list
    msg = MIMEMultipart()
    msg['From'] = send_from
    msg['To'] = send_to
    msg['Date'] = formatdate(localtime = True)
    if not path_found:
        msg['Subject'] = "###CUMULATIVE STATS NOT UPDATED - Script1 stats summary"
    else:
        msg['Subject'] = "Script1 stats summary"
    msg.attach(MIMEText(text))

    part = MIMEBase('application', "octet-stream")
    part.set_payload(open(summary_file, "rb").read())
    encoders.encode_base64(part)
    part.add_header('Content-Disposition', 'attachment; filename="summary_file.xlsx"')
    msg.attach(part)

    #context = ssl.SSLContext(ssl.PROTOCOL_SSLv3)
    #SSL connection only working on Python 3+
    smtp = smtplib.SMTP('10.10.8.12')

    smtp.send_message(msg)
    #smtp.sendmail(send_from, send_to, msg.as_string())
    smtp.quit()

if options.email:
    send_email(email_list)


print ("\n\nruntime: %s:  \n" % runtime)






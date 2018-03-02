import json
import xlsxwriter
import os
import shutil
import pandas as pd

root_dir = str(os.getcwd())

def step1_stats_out(master_stat_summary):

    with open(master_stat_summary) as infile:
        master_stat_list_dict = json.load(infile)

    st = "test"
    summary_file = root_dir+ '/stat_alignment_summary_' + st + '.xlsx'
    workbook = xlsxwriter.Workbook(summary_file)
    worksheet = workbook.add_worksheet()
    row = 0
    col = 0
    top_row_header = ["time_stamp", "sample_name", "self.species", "reference_sequence_name", "R1size", "R2size", "allbam_mapped_reads", "genome_coverage", "ave_coverage", "ave_read_length", "unmapped_reads", "unmapped_assembled_contigs", "good_snp_count", "mlst_type", "octalcode", "sbcode", "hexadecimal_code", "binarycode"]
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
        summary_cumulative_file_temp = copy_to + '/stat_alignment_culmulative_summary-' + st + '-temp.xlsx'
        temp_folder = copy_to + '/temp'

    frames = []
    for stat_summary_dict in master_stat_list_dict:
        df_stat_summary = pd.DataFrame.from_dict(stat_summary_dict, orient='index') #convert stat_summary to df
        frames.append(df_stat_summary) #frames to concatenate
        col = 0
        row += 1
        #stat summary to be attached in email and left in root directory
        for v in stat_summary_dict.values():
            worksheet.write(row, col, v)
            col += 1

    workbook.close()

    if path_found:
        try:
            open_check = open(summary_cumulative_file, 'a') #'a' is very important, 'w' will leave you with an empty file
            open_check.close()
            df_all=pd.read_excel(summary_cumulative_file)
            df_all_trans = df_all.T #indexed on column headers
            # save back the old and remake the working stats file
            shutil.move(summary_cumulative_file, '{}' .format(temp_folder + '/stat_backup' + st + '.xlsx'))
            sorter = list(df_all_trans.index) #list of original column order
            frames.insert(0, df_all_trans) #put as first item in list
            df_concat = pd.concat(frames, axis=1) #cat frames
            df_sorted = df_concat.loc[sorter] #sort based on sorter order
            df_sorted.T.to_excel(summary_cumulative_file, index=False) #transpose before writing to excel, numerical index not needed
        except BlockingIOError:
            sorter = list(df_stat_summary.index) #list of original column order
            df_concat = pd.concat(frames, axis=1) #cat frames
            df_sorted = df_concat.loc[sorter] #sort based on sorter order
            df_sorted.T.to_excel(summary_cumulative_file_temp, index=False)
        except OSError:
            sorter = list(df_stat_summary.index) #list of original column order
            df_concat = pd.concat(frames, axis=1) #cat frames
            df_sorted = df_concat.loc[sorter] #sort based on sorter order
            df_sorted.T.to_excel(summary_cumulative_file_temp, index=False)
    else:
        print("Path to cumulative stat summary file not found")






def iterating_dirs(run_list):
             frames = []
             stat_summary

                 df_stat_summary = pd.DataFrame.from_dict(stat_summary, orient='index') #convert stat_summary to df
                 frames.append(df_stat_summary) #frames to concatenate
                     col = 0
                     row += 1
                     #run stats
                     for v in stat_summary.values():
                         worksheet.write(row, col, v) #stat summary to be attached in email and left in root directory
                         col += 1
             workbook.close()
             if not quiet_call and path_found:
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


 if args.email == "none":
     print ("\n\temail not sent")
 elif args.email:
     send_email()
     print ("\n\temail sent to: %s" % email_list)
 else:
     print ("\n\temail not sent")

 if args.upload:
     print ("Uploading Samples...")

     #upload to bioinfoVCF
     src = root_dir
     dst = bioinfoVCF + "/" + os.path.basename(os.path.normpath(root_dir))
     print ("\n\t%s is copying to %s" % (src, dst))
     os.makedirs(dst, exist_ok=True)
     copy_tree(src, dst, preserve_mode=0, preserve_times=0)


import os
import shutil
import git
import glob
import re
import xlsxwriter
import xlrd

class Update_Directory:
    def __init__(self, dependents_dir):
        home = os.path.expanduser("~")

        if os.path.isdir("/bioinfo11/TStuber/Results"):
            #server
            storage_path = "/bioinfo11/TStuber/Results"
            computer_path = "/home/shared"
        elif os.path.isdir("/Volumes/root/TStuber/Results"):
            #mac
            storage_path = "/Volumes/root/TStuber/Results"
            computer_path = "/Users/Shared"

        #check if bioinfo is attached
        if os.path.isdir(storage_path): 
            self.upload_to = storage_path
            self.remote = storage_path + dependents_dir
            self.path = storage_path + dependents_dir
            #make copy of local dependents
            if os.path.isdir(computer_path):
                dep_path = computer_path
                dir_split = dependents_dir.split('/')[1:]
                #build path to file if needed
                for i in dir_split:
                    dep_path += '/' + i
                    if not os.path.exists(dep_path):
                        os.makedirs(dep_path)
                local = computer_path + dependents_dir
                if os.path.isdir(local):
                    try:
                        #remove files
                        shutil.rmtree(local)
                        #copy files from remote to local
                        shutil.copytree(self.remote, local)
                    except:
                        pass
        #if no storage path (bioinfo) check if dependents have been copied local to share
        elif os.path.isdir(computer_path + dependents_dir):
            self.upload_to = "not_found"
            self.remote = "not_found"
            self.path = computer_path + dependents_dir
        #if not in shared folder then check home directory previously downloaded from github
        elif os.path.isdir(home + "/dependencies" + dependents_dir):
            self.upload_to = "not_found"
            self.remote = "not_found"
            self.path = home + "/dependencies" + dependents_dir # sets dependencies directory to home directory
        else:
            os.makedirs(home + "/dependencies")
            print("\n\nDOWNLOADING DEPENDENCIES FROM GITHUB... ***\n\n")
            git.Repo.clone_from("https://github.com/USDA-VS/dependencies.git", home + "/dependencies")
            self.upload_to = "not_found"
            self.remote = "not_found"
            self.path = home + "/dependencies" + dependents_dir

def Bruc_Private_Codes(genotypingcodes):
    found = False
    if os.path.isfile("/Volumes/MB/Brucella/Brucella Logsheets/ALL_WGS.xlsx"):
        private_location = "/Volumes/MB/Brucella/Brucella Logsheets/ALL_WGS.xlsx"
        print("Copying single column genotype codes to:  %s" % genotypingcodes)
        found = True

    elif os.path.isfile("/fdrive/Brucella/Brucella Logsheets/ALL_WGS.xlsx"):
        private_location = "/fdrive/Brucella/Brucella Logsheets/ALL_WGS.xlsx"
        print("Copying single column genotype codes to:  %s" % genotypingcodes)
        found = True

    else:
        print("Path to Brucella genotyping codes not found")

    if found:
        wb_out = xlsxwriter.Workbook(genotypingcodes)
        ws_out = wb_out.add_worksheet()

        wb_in = xlrd.open_workbook(private_location)

        row = 0
        col = 0

        sheet_in = wb_in.sheet_by_index(1)
        for row_data in sheet_in.col(32):
            row_data = row_data.value
            row_data = re.sub("/", "_", row_data)
            row_data = re.sub("\.", "_", row_data)
            row_data = re.sub("\*", "_", row_data)
            row_data = re.sub("\?", "_", row_data)
            row_data = re.sub("\(", "_", row_data)
            row_data = re.sub("\)", "_", row_data)
            row_data = re.sub("\[", "_", row_data)
            row_data = re.sub("\]", "_", row_data)
            row_data = re.sub(" ", "_", row_data)
            row_data = re.sub("{", "_", row_data)
            row_data = re.sub("}", "_", row_data)
            row_data = re.sub("\'", "_", row_data)
            row_data = re.sub("-_", "_", row_data)
            row_data = re.sub("_-", "_", row_data)
            row_data = re.sub("--", "_", row_data)
            row_data = re.sub("_$", "", row_data)
            row_data = re.sub("-$", "", row_data)
            row_data = re.sub("\'", "", row_data)
            row_data = str(row_data)

            ws_out.write(row, col, row_data)
            row += 1

        wb_out.close()

def Get_Filters(self, excelinfile, filter_files):
    for i in glob.glob(filter_files + "/*"):
        os.remove(i)

    wb = xlrd.open_workbook(excelinfile)
    sheets = wb.sheet_names()
    for sheet in sheets:
        ws = wb.sheet_by_name(sheet)

        myrange = lambda start, end: range(start, end+1)

        for colnum in range(ws.ncols): # for each column in worksheet
            file_out = filter_files + "/" + ws.col_values(colnum)[0] + ".txt" # column header naming file
            write_out = open (file_out, 'at')
            mylist = ws.col_values(colnum)[1:] # list of each field in column, minus the header
            mylist = [x for x in mylist if x] # remove blank cells
            for value in mylist:
                value = str(value)
                value = value.replace(sheet + "-", '')
                if "-" not in value:
                    value=int(float(value)) # change str to float to int
                    print (sheet + "-" + str(value), file=write_out)
                elif "-" in value:
                    value = value.split("-")
                    for i in range(int(value[0]), int(value[1]) + 1 ):
                        print (sheet + "-" + str(i), file=write_out)
    write_out.close()
import os
import shutil
import sys
import glob
import git
import xlrd

class Update_Directory:
    def __init__(self, dependents_dir):
        home = os.path.expanduser("~")
        storage_path = ""

        if os.path.isdir("/bioinfo11/TStuber/Results"):
            #server
            storage_path = "/bioinfo11/TStuber/Results"
            computer_path = "/home/shared"
        elif os.path.isdir("/Volumes/root/TStuber/Results"):
            #mac
            storage_path = "/Volumes/root/TStuber/Results"
            computer_path = "/Users/Shared"
        elif os.path.isdir("/home/shared"):
            computer_path = "/home/shared"
        elif os.path.isdir("/Users/Shared"):
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


class Step1:

    '''
        Internally for each species step1_dependents is made to an object.  This returns 3 attributes: upload_to, remote and path 

        from vAttributes import Step1
        Instantiation: x = Step1("ab1"), x.reference
            provides 7 attributes

    '''
    def __init__(self, species):
        if species == "ab1":
            self.dependents_dir="/brucella/abortus1/script_dependents/script1"
            step1_dependents = Update_Directory(self.dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/NC_00693c.fasta"
            print("Reference used: %s" % self.reference)
            self.hqs = step1_dependents.path + "/NC_00693cHighestQualitySNPs.vcf"
            self.gbk_file = step1_dependents.path + "/NC_006932-NC_006933.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
        
        elif species == "ab3":
            self.dependents_dir="/brucella/abortus3/script_dependents/script1"
            step1_dependents = Update_Directory(self.dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/CP007682-7683c.fasta"
            print("Reference used: %s" % self.reference)
            self.hqs = step1_dependents.path + "/CP007682-7683cHighestQualitySNPs.vcf"
            self.gbk_file = step1_dependents.path + "/CP007682-CP007683.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
            
        elif species == "canis":
            self.dependents_dir="/brucella/canis/script_dependents/script1"
            step1_dependents = Update_Directory(self.dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/BcanisATCC23365.fasta"
            print("Reference used: %s" % self.reference)
            self.hqs = step1_dependents.path + "/canisHighestQualitySNPs.vcf"
            self.gbk_file = step1_dependents.path + "/NC_010103-NC_010104.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
        
        elif species == "ceti1":
            self.dependents_dir="/brucella/ceti1/script_dependents/script1"
            step1_dependents = Update_Directory(self.dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/Bceti1Cudo.fasta"
            print("Reference used: %s" % self.reference)
            self.hqs = step1_dependents.path + "/ceti1HighestQualitySNPs.vcf"
            self.gbk_file = None #script_dependents + "/no.gff"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path

        elif species == "ceti2":
            self.dependents_dir="/brucella/ceti2/script_dependents/script1"
            step1_dependents = Update_Directory(self.dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/Bceti2-TE10759.fasta"
            print("Reference used: %s" % self.reference)
            self.hqs = step1_dependents.path + "/ceti2HighestQualitySNPs.vcf"
            self.gbk_file = step1_dependents.path + "/NC_022905-NC_022906.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path

        elif species == "mel1":
            self.dependents_dir="/brucella/melitensis-bv1/script_dependents/script1"
            step1_dependents = Update_Directory(self.dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/mel-bv1-NC003317.fasta"
            print("Reference used: %s" % self.reference)
            self.hqs = step1_dependents.path + "/mel-bv1-NC003317-highqualitysnps.vcf"
            self.gbk_file = step1_dependents.path + "/NC_003317-NC_003318.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
        
        elif species == "mel1b":
            self.dependents_dir="/brucella/melitensis-bv1b/script_dependents/script1"
            step1_dependents = Update_Directory(self.dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/mel-bv1b-CP018508.fasta"
            print("Reference used: %s" % self.reference)
            self.hqs = step1_dependents.path + "/mel-bv1b-CP018508-highqualitysnps.vcf"
            self.gbk_file = step1_dependents.path + "/mel-bv1b-CP018508.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
        
        elif species == "mel2":
            self.dependents_dir="/brucella/melitensis-bv2/script_dependents/script1"
            step1_dependents = Update_Directory(self.dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/mel-bv2-NC012441.fasta"
            print("Reference used: %s" % self.reference)
            self.hqs = step1_dependents.path + "/mel-bv2-NC012441-highqualitysnps.vcf"
            self.gbk_file = step1_dependents.path + "/NC_012441-NC_012442.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
            
        elif species == "mel3":
            self.dependents_dir="/brucella/melitensis-bv3/script_dependents/script1"
            step1_dependents = Update_Directory(self.dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/mel-bv3-NZCP007760.fasta"
            print("Reference used: %s" % self.reference)
            self.hqs = step1_dependents.path + "/mel-bv3-NZCP007760-highqualitysnps.vcf"
            self.gbk_file = step1_dependents.path + "/NZ_CP007760-NZ_CP007761.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
        
        elif species == "suis1":
            self.dependents_dir="/brucella/suis1/script_dependents/script1"
            step1_dependents = Update_Directory(self.dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/NC_017251-NC_017250.fasta"
            print("Reference used: %s" % self.reference)
            self.hqs = step1_dependents.path + "/B00-0468-highqualitysnps.vcf"
            self.gbk_file = step1_dependents.path + "/NC_017251-NC_017250.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
            
        elif species == "suis3":
            self.dependents_dir="/brucella/suis3/script_dependents/script1"
            step1_dependents = Update_Directory(self.dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/NZ_CP007719-NZ_CP007718.fasta"
            print("Reference used: %s" % self.reference)
            self.hqs = step1_dependents.path + "/highqualitysnps.vcf"
            self.gbk_file = None #script_dependents + "/no.gff"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
            
        elif species == "suis4":
            self.dependents_dir="/brucella/suis4/script_dependents/script1"
            step1_dependents = Update_Directory(self.dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/B-REF-BS4-40.fasta"
            print("Reference used: %s" % self.reference)
            self.hqs = step1_dependents.path + "/suis4HighestQualitySNPs.vcf"
            self.gbk_file = None #script_dependents + "/no.gff"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
       
        elif species == "ovis":
            self.dependents_dir="/brucella/ovis/script_dependents/script1"
            step1_dependents = Update_Directory(self.dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/BovisATCC25840.fasta"
            print("Reference used: %s" % self.reference)
            self.hqs = step1_dependents.path + "/BovisATCC25840HighestQualitySNPs.vcf"
            self.gbk_file = step1_dependents.path + "/NC_009505-NC_009504.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
            
        elif species == "neo":
            self.dependents_dir="/brucella/neotomae/script_dependents/script1"
            step1_dependents = Update_Directory(self.dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/KN046827.fasta"
            print("Reference used: %s" % self.reference)
            self.hqs = step1_dependents.path + "/ERR1845155-highqualitysnps.vcf"
            self.gbk_file = step1_dependents.path + "/KN046827.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
            
        elif species == "af":
            self.dependents_dir="/mycobacterium/tbc/af2122/script_dependents/script1"
            step1_dependents = Update_Directory(self.dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/spoligotype_db.txt"
            self.reference = step1_dependents.path + "/NC_002945v4.fasta"
            print("Reference used: %s" % self.reference)
            self.hqs = step1_dependents.path + "/highqualitysnps.vcf"
            self.gbk_file = step1_dependents.path + "/NC_002945v4.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
        
        elif species == "h37":
            self.dependents_dir="/mycobacterium/tbc/h37/script_dependents/script1"
            step1_dependents = Update_Directory(self.dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/spoligotype_db.txt"
            self.reference = step1_dependents.path + "/NC000962.fasta"
            print("Reference used: %s" % self.reference)
            self.hqs = step1_dependents.path + "/15-3162-highqualitysnps.vcf"
            self.gbk_file = step1_dependents.path + "/NC_000962.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
            
        elif species == "para":
            self.dependents_dir="/mycobacterium/avium_complex/para_cattle-bison/script_dependents/script1"
            step1_dependents = Update_Directory(self.dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/NC_002944.fasta"
            print("Reference used: %s" % self.reference)
            self.hqs = step1_dependents.path + "/HQ-NC002944.vcf"
            self.gbk_file = step1_dependents.path + "/NC_002944.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
            
class Step2:
    def __init__(self, species):
                            
        if species == "suis1":
            qual_gatk_threshold = 300
            N_gatk_threshold = 350
            dependents_dir="/brucella/suis1/script_dependents/script2"
            step2_dependents = Update_Directory(dependents_dir)
           
            #step2_dependents.upload_to --> '/Volumes/root/TStuber/Results'
            genotypingcodes = step2_dependents.upload_to + "/brucella/genotyping_codes.xlsx" # this may not be available if there is no access to f drive.  
            Bruc_Private_Codes(genotypingcodes) # if f drive then upload fixed column 32 to bioinfo
            
            gbk_file = step2_dependents.path + "/NC_017251-NC_017250.gbk"
            definingSNPs = step2_dependents.path + "/DefiningSNPsGroupDesignations_python.xlsx"
            remove_from_analysis = step2_dependents.path + "/RemoveFromAnalysis.xlsx"
            bioinfoVCF = step2_dependents.upload_to + "/brucella/suis1/vcfs"
            excelinfile = step2_dependents.path + "/Filtered_Regions_python.xlsx"

            #Get_Filters(excelinfile) #***FUNCTION CALL

        elif species == "suis3":
            qual_gatk_threshold = 300
            N_gatk_threshold = 350
            dependents_dir="/brucella/suis3/script_dependents/script2"
            step2_dependents = Update_Directory(dependents_dir)

            genotypingcodes = step2_dependents.upload_to + "/brucella/genotyping_codes.xlsx"
            Bruc_Private_Codes(genotypingcodes)
            
            gbk_file = step2_dependents.path + "/NZ_CP007719-NZ_CP007718.gbk"
            definingSNPs = step2_dependents.path + "/DefiningSNPsGroupDesignations_python.xlsx"
            remove_from_analysis = step2_dependents.path + "/RemoveFromAnalysis.xlsx"
            bioinfoVCF = step2_dependents.upload_to + "/brucella/suis3/vcfs"
            excelinfile = step2_dependents.path + "/Filtered_Regions_python.xlsx"

        elif species == "suis4":
            qual_gatk_threshold = 300
            N_gatk_threshold = 350
            dependents_dir="/brucella/suis4/script_dependents/script2"
            step2_dependents = Update_Directory(dependents_dir)

            genotypingcodes = step2_dependents.upload_to + "/brucella/genotyping_codes.xlsx"
            Bruc_Private_Codes(genotypingcodes)
            
            #gbk_file = step2_dependents.path + ""
            definingSNPs = step2_dependents.path + "/DefiningSNPsGroupDesignations_python.xlsx"
            remove_from_analysis = step2_dependents.path + "/RemoveFromAnalysis.xlsx"
            bioinfoVCF = step2_dependents.upload_to + "/brucella/suis4/vcfs"
            excelinfile = step2_dependents.path + "/Filtered_Regions_python.xlsx"
 
        elif species == "ab1":
            qual_gatk_threshold = 300
            N_gatk_threshold = 350
            dependents_dir="/brucella/abortus1/script_dependents/script2"
            step2_dependents = Update_Directory(dependents_dir)

            genotypingcodes = step2_dependents.upload_to + "/brucella/genotyping_codes.xlsx"
            Bruc_Private_Codes(genotypingcodes)
            
            gbk_file = step2_dependents.path + "/NC_006932-NC_006933.gbk"
            definingSNPs = step2_dependents.path + "/DefiningSNPsGroupDesignations_python.xlsx"
            remove_from_analysis = step2_dependents.path + "/RemoveFromAnalysis.xlsx"
            bioinfoVCF = step2_dependents.upload_to + "/brucella/abortus1/vcfs"
            excelinfile = step2_dependents.path + "/Filtered_Regions_python.xlsx"

        elif species == "ab3":
            qual_gatk_threshold = 300
            N_gatk_threshold = 350
            dependents_dir="/brucella/abortus3/script_dependents/script2"
            step2_dependents = Update_Directory(dependents_dir)

            genotypingcodes = step2_dependents.upload_to + "/brucella/genotyping_codes.xlsx"
            Bruc_Private_Codes(genotypingcodes)
            
            gbk_file = step2_dependents.path + "/CP007682-CP007683.gbk"
            definingSNPs = step2_dependents.path + "/DefiningSNPsGroupDesignations_python.xlsx"
            remove_from_analysis = step2_dependents.path + "/RemoveFromAnalysis.xlsx"
            bioinfoVCF = step2_dependents.upload_to + "/brucella/abortus3/vcfs"
            excelinfile = step2_dependents.path + "/Filtered_Regions_python.xlsx"
            print(excelinfile)
 
        elif species == "mel1":
            qual_gatk_threshold = 300
            N_gatk_threshold = 350
            dependents_dir="/brucella/melitensis-bv1/script_dependents/script2"
            step2_dependents = Update_Directory(dependents_dir)

            genotypingcodes = step2_dependents.upload_to + "/brucella/genotyping_codes.xlsx"
            Bruc_Private_Codes(genotypingcodes)
            
            gbk_file = step2_dependents.path + "/NC_003317-NC_003318.gbk"
            definingSNPs = step2_dependents.path + "/DefiningSNPsGroupDesignations_python.xlsx"
            remove_from_analysis = step2_dependents.path + "/RemoveFromAnalysis.xlsx"
            bioinfoVCF = step2_dependents.upload_to + "/brucella/melitensis-bv1/vcfs"
            excelinfile = step2_dependents.path + "/Filtered_Regions_python.xlsx"
            print(excelinfile)

        elif species == "mel1b":
            qual_gatk_threshold = 300
            N_gatk_threshold = 350
            dependents_dir="/brucella/melitensis-bv1b/script_dependents/script2"
            step2_dependents = Update_Directory(dependents_dir)

            genotypingcodes = step2_dependents.upload_to + "/brucella/genotyping_codes.xlsx"
            Bruc_Private_Codes(genotypingcodes)
            
            gbk_file = step2_dependents.path + "/mel-bv1b-CP018508.gbk"
            definingSNPs = step2_dependents.path + "/DefiningSNPsGroupDesignations.xlsx"
            remove_from_analysis = step2_dependents.path + "/RemoveFromAnalysis.xlsx"
            bioinfoVCF = step2_dependents.upload_to + "/brucella/melitensis-bv1b/vcfs"
            excelinfile = step2_dependents.path + "/Filtered_Regions.xlsx"
            print(excelinfile)

        elif species == "mel2":
            qual_gatk_threshold = 300
            N_gatk_threshold = 350
            dependents_dir="/brucella/melitensis-bv2/script_dependents/script2"
            step2_dependents = Update_Directory(dependents_dir)

            genotypingcodes = step2_dependents.upload_to + "/brucella/genotyping_codes.xlsx"
            Bruc_Private_Codes(genotypingcodes)
            
            gbk_file = step2_dependents_path + "/NC_012441-NC_012442.gbk"
            definingSNPs = step2_dependents.path + "/DefiningSNPsGroupDesignations_python.xlsx"
            remove_from_analysis = step2_dependents.path + "/RemoveFromAnalysis.xlsx"
            bioinfoVCF = step2_dependents.upload_to + "/brucella/melitensis-bv2/vcfs"
            excelinfile = step2_dependents.path + "/Filtered_Regions_python.xlsx"
            print(excelinfile)

        elif species == "mel3":
            qual_gatk_threshold = 300
            N_gatk_threshold = 350
            dependents_dir="/brucella/melitensis-bv3/script_dependents/script2"
            step2_dependents = Update_Directory(dependents_dir)

            genotypingcodes = step2_dependents.upload_to + "/brucella/genotyping_codes.xlsx"
            Bruc_Private_Codes(genotypingcodes)
            
            gbk_file = step2_dependents.path + "/NZ_CP007760-NZ_CP007761.gbk"
            definingSNPs = step2_dependents.path + "/DefiningSNPsGroupDesignations_python.xlsx"
            remove_from_analysis = step2_dependents.path + "/RemoveFromAnalysis.xlsx"
            bioinfoVCF = step2_dependents.upload_to + "/brucella/melitensis-bv3/vcfs"
            excelinfile = step2_dependents.path + "/Filtered_Regions_python.xlsx"
            print(excelinfile)

        elif species == "canis":
            qual_gatk_threshold = 300
            N_gatk_threshold = 350
            dependents_dir="/brucella/canis/script_dependents/script2"
            step2_dependents = Update_Directory(dependents_dir)

            genotypingcodes = step2_dependents.upload_to + "/brucella/genotyping_codes.xlsx"
            Bruc_Private_Codes(genotypingcodes)
            
            gbk_file = step2_dependents.path + "/NC_010103-NC_010104.gbk"
            definingSNPs = step2_dependents.path + "/DefiningSNPsGroupDesignations_python.xlsx"
            remove_from_analysis = step2_dependents.path + "/RemoveFromAnalysis.xlsx"
            bioinfoVCF = step2_dependents.upload_to + "/brucella/canis/vcfs"
            excelinfile = step2_dependents.path + "/Filtered_Regions_python.xlsx"
            print(excelinfile)

        elif species == "ceti1":
            qual_gatk_threshold = 300
            N_gatk_threshold = 350
            dependents_dir="/brucella/ceti1/script_dependents/script2"
            step2_dependents = Update_Directory(dependents_dir)

            genotypingcodes = step2_dependents.upload_to + "/brucella/genotyping_codes.xlsx"
            Bruc_Private_Codes(genotypingcodes)
            
            #gbk_file = step1_dependents_path + ""
            definingSNPs = step2_dependents.path + "/DefiningSNPsGroupDesignations_python.xlsx"
            remove_from_analysis = step2_dependents.path + "/RemoveFromAnalysis.xlsx"
            bioinfoVCF = step2_dependents.upload_to + "/brucella/ceti1/vcfs"
            excelinfile = step2_dependents.path + "/Filtered_Regions.xlsx"
            print(excelinfile)

        elif species == "ceti2":
            qual_gatk_threshold = 300
            N_gatk_threshold = 350
            dependents_dir="/brucella/ceti2/script_dependents/script2"
            step2_dependents = Update_Directory(dependents_dir)

            genotypingcodes = step2_dependents.upload_to + "/brucella/genotyping_codes.xlsx"
            Bruc_Private_Codes(genotypingcodes)
            
            gbk_file = step2_dependents.path + "/NC_022905-NC_022906.gbk"
            definingSNPs = step2_dependents.path + "/DefiningSNPsGroupDesignations_python.xlsx"
            remove_from_analysis = step2_dependents.path + "/RemoveFromAnalysis.xlsx"
            bioinfoVCF = step2_dependents.upload_to + "/brucella/ceti2/vcfs"
            excelinfile = step2_dependents.path + "/Filtered_Regions.xlsx"
            print(excelinfile)
                
        elif species == "ovis":
            qual_gatk_threshold = 300
            N_gatk_threshold = 350
            dependents_dir="/brucella/ovis/script_dependents/script2"
            step2_dependents = Update_Directory(dependents_dir)

            genotypingcodes = step2_dependents.upload_to + "/brucella/genotyping_codes.xlsx"
            Bruc_Private_Codes(genotypingcodes)
            
            gbk_file = step2_dependents.path + "/NC_009505-NC_009504.gbk"
            definingSNPs = step2_dependents.path + "/DefiningSNPsGroupDesignations.xlsx"
            remove_from_analysis = step2_dependents.path + "/RemoveFromAnalysis.xlsx"
            bioinfoVCF = step2_dependents.upload_to + "/brucella/ovis/vcfs"
            excelinfile = step2_dependents.path + "/Filtered_Regions.xlsx"
            print(excelinfile)
                
        elif species == "neo":
            qual_gatk_threshold = 300
            N_gatk_threshold = 350
            dependents_dir="/brucella/neotomae/script_dependents/script2"
            step2_dependents = Update_Directory(dependents_dir)

            genotypingcodes = step2_dependents.upload_to + "/brucella/genotyping_codes.xlsx"
            Bruc_Private_Codes(genotypingcodes)
            
            gbk_file = step2_dependents.path + "/KN046827.gbk"
            definingSNPs = step2_dependents.path + "/DefiningSNPsGroupDesignations.xlsx"
            remove_from_analysis = step2_dependents.path + "/RemoveFromAnalysis.xlsx"
            bioinfoVCF = step2_dependents.upload_to + "/brucella/neotomae/vcfs"
            excelinfile = step2_dependents.path + "/Filtered_Regions.xlsx"
            print(excelinfile)

        elif species == "af":
            qual_gatk_threshold = 150
            N_gatk_threshold = 200
            dependents_dir="/mycobacterium/tbc/af2122/script_dependents/script2"
            step2_dependents = Update_Directory(dependents_dir)

            genotypingcodes = step2_dependents.upload_to + "/mycobacterium/genotyping_codes.xlsx"

            gbk_file = step2_dependents.path + "/NC_002945v4.gbk"
            definingSNPs = step2_dependents.path + "/DefiningSNPsGroupDesignations.xlsx"
            remove_from_analysis = step2_dependents.path + "/RemoveFromAnalysis.xlsx"
            bioinfoVCF = step2_dependents.upload_to + "/mycobacterium/tbc/af2122/script2"
            excelinfile = step2_dependents.path + "/Filtered_Regions.xlsx"

        elif species == "h37":
            qual_gatk_threshold = 150
            N_gatk_threshold = 200
            dependents_dir="/mycobacterium/tbc/h37/script_dependents/script2"
            step2_dependents = Update_Directory(dependents_dir)

            genotypingcodes = step2_dependents.upload_to + "/mycobacterium/genotyping_codes.xlsx"
            
            gbk_file = step2_dependents.path + "/NC_000962.gbk"
            definingSNPs = step2_dependents.path + "/DefiningSNPsGroupDesignations_python.xlsx"
            remove_from_analysis = step2_dependents.path + "/RemoveFromAnalysis.xlsx"
            bioinfoVCF = step2_dependents.upload_to + "/mycobacterium/tbc/h37/script2"
            excelinfile = step2_dependents.path + "/Filtered_Regions_python.xlsx"

        elif species == "para":
            qual_gatk_threshold = 150
            N_gatk_threshold = 200
            dependents_dir="/mycobacterium/avium_complex/para_cattle-bison/script_dependents/script2"
            step2_dependents = Update_Directory(dependents_dir)

            genotypingcodes = step2_dependents.upload_to + "/mycobacterium/genotyping_codes.xlsx"
            
            gbk_file = step2_dependents.path + "/NC_002944.gbk"
            definingSNPs = step2_dependents.path + "/DefiningSNPsGroupDesignations.xlsx"
            remove_from_analysis = step2_dependents.path + "/RemoveFromAnalysis.xlsx"
            bioinfoVCF = step2_dependents.upload_to + "/mycobacterium/avium_complex/para_cattle-bison/vcfs"
            excelinfile = step2_dependents.path + "/Filtered_Regions.xlsx"

        else:
            print ("\n#####EXIT AT SETTING OPTIONS, Check that a \"-s\" species was provided\n")
            sys.exit(0)
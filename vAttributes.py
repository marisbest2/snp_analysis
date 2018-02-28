import os
import shutil
import sys

# vSNP classes
from vFunctions import Update_Directory
from vFunctions import Bruc_Private_Codes
from vFunctions import Get_Filters

class Step1:

    '''
        Internally for each species step1_dependents is made to an object.  This returns 3 attributes: upload_to, remote and path 

        from vAttributes import Step1
        Instantiation: x = Step1("ab1"), x.reference
            provides 10 attributes

'''
    def __init__(self, species):
        if species == "ab1":
            self.dependents_dir="/brucella/abortus1/script_dependents/script1"
            self.step1_dependents = Update_Directory(dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/NC_00693c.fasta"
            print("Reference being used: %s" % reference)
            self.hqs = step1_dependents.path + "/NC_00693cHighestQualitySNPs.vcf"
            self.gbk_file = step1_dependents.path + "/NC_006932-NC_006933.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
        
        if species == "ab3":
            self.dependents_dir="/brucella/abortus3/script_dependents/script1"
            self.step1_dependents = Update_Directory(dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/CP007682-7683c.fasta"
            print("Reference being used: %s" % reference)
            self.hqs = step1_dependents.path + "/CP007682-7683cHighestQualitySNPs.vcf"
            self.gbk_file = step1_dependents.path + "/CP007682-CP007683.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
            
        if species == "canis":
            self.dependents_dir="/brucella/canis/script_dependents/script1"
            self.step1_dependents = Update_Directory(dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/BcanisATCC23365.fasta"
            print("Reference being used: %s" % reference)
            self.hqs = step1_dependents.path + "/canisHighestQualitySNPs.vcf"
            self.gbk_file = step1_dependents.path + "/NC_010103-NC_010104.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
        
        if species == "ceti1":
            self.dependents_dir="/brucella/ceti1/script_dependents/script1"
            self.step1_dependents = Update_Directory(dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/Bceti1Cudo.fasta"
            print("Reference being used: %s" % reference)
            self.hqs = step1_dependents.path + "/ceti1HighestQualitySNPs.vcf"
            self.gbk_file = None #script_dependents + "/no.gff"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path

        if species == "ceti2":
            self.dependents_dir="/brucella/ceti2/script_dependents/script1"
            self.step1_dependents = Update_Directory(dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/Bceti2-TE10759.fasta"
            print("Reference being used: %s" % reference)
            self.hqs = step1_dependents.path + "/ceti2HighestQualitySNPs.vcf"
            self.gbk_file = step1_dependents.path + "/NC_022905-NC_022906.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path

        if species == "mel1":
            self.dependents_dir="/brucella/melitensis-bv1/script_dependents/script1"
            self.step1_dependents = Update_Directory(dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/mel-bv1-NC003317.fasta"
            print("Reference being used: %s" % reference)
            self.hqs = step1_dependents.path + "/mel-bv1-NC003317-highqualitysnps.vcf"
            self.gbk_file = step1_dependents.path + "/NC_003317-NC_003318.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
        
        if species == "mel1b":
            self.dependents_dir="/brucella/melitensis-bv1b/script_dependents/script1"
            self.step1_dependents = Update_Directory(dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/mel-bv1b-CP018508.fasta"
            print("Reference being used: %s" % reference)
            self.hqs = step1_dependents.path + "/mel-bv1b-CP018508-highqualitysnps.vcf"
            self.gbk_file = step1_dependents.path + "/mel-bv1b-CP018508.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
        
        if species == "mel2":
            self.dependents_dir="/brucella/melitensis-bv2/script_dependents/script1"
            self.step1_dependents = Update_Directory(dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/mel-bv2-NC012441.fasta"
            print("Reference being used: %s" % reference)
            self.hqs = step1_dependents.path + "/mel-bv2-NC012441-highqualitysnps.vcf"
            self.gbk_file = step1_dependents.path + "/NC_012441-NC_012442.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
            
        if species == "mel3":
            self.dependents_dir="/brucella/melitensis-bv3/script_dependents/script1"
            self.step1_dependents = Update_Directory(dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/mel-bv3-NZCP007760.fasta"
            print("Reference being used: %s" % reference)
            self.hqs = step1_dependents.path + "/mel-bv3-NZCP007760-highqualitysnps.vcf"
            self.gbk_file = step1_dependents.path + "/NZ_CP007760-NZ_CP007761.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
        
        if species == "suis1":
            self.dependents_dir="/brucella/suis1/script_dependents/script1"
            self.step1_dependents = Update_Directory(dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/NC_017251-NC_017250.fasta"
            print("Reference being used: %s" % reference)
            self.hqs = step1_dependents.path + "/B00-0468-highqualitysnps.vcf"
            self.gbk_file = step1_dependents.path + "/NC_017251-NC_017250.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
            
        if species == "suis3":
            self.dependents_dir="/brucella/suis3/script_dependents/script1"
            self.step1_dependents = Update_Directory(dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/NZ_CP007719-NZ_CP007718.fasta"
            print("Reference being used: %s" % reference)
            self.hqs = step1_dependents.path + "/highqualitysnps.vcf"
            self.gbk_file = None #script_dependents + "/no.gff"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
            
        if species == "suis4":
            self.dependents_dir="/brucella/suis4/script_dependents/script1"
            self.step1_dependents = Update_Directory(dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/B-REF-BS4-40.fasta"
            print("Reference being used: %s" % reference)
            self.hqs = step1_dependents.path + "/suis4HighestQualitySNPs.vcf"
            self.gbk_file = None #script_dependents + "/no.gff"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
       
        if species == "ovis":
            self.dependents_dir="/brucella/ovis/script_dependents/script1"
            self.step1_dependents = Update_Directory(dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/BovisATCC25840.fasta"
            print("Reference being used: %s" % reference)
            self.hqs = step1_dependents.path + "/BovisATCC25840HighestQualitySNPs.vcf"
            self.gbk_file = step1_dependents.path + "/NC_009505-NC_009504.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
            
        if species == "neo":
            self.dependents_dir="/brucella/neotomae/script_dependents/script1"
            self.step1_dependents = Update_Directory(dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/KN046827.fasta"
            print("Reference being used: %s" % reference)
            self.hqs = step1_dependents.path + "/ERR1845155-highqualitysnps.vcf"
            self.gbk_file = step1_dependents.path + "/KN046827.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
            
        if species == "af":
            self.dependents_dir="/mycobacterium/tbc/af2122/script_dependents/script1"
            self.step1_dependents = Update_Directory(dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/spoligotype_db.txt"
            self.reference = step1_dependents.path + "/NC_002945v4.fasta"
            print("Reference being used: %s" % reference)
            self.hqs = step1_dependents.path + "/highqualitysnps.vcf"
            self.gbk_file = step1_dependents.path + "/NC_002945v4.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
        
        if species == "h37":
            self.dependents_dir="/mycobacterium/tbc/h37/script_dependents/script1"
            self.step1_dependents = Update_Directory(dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/spoligotype_db.txt"
            self.reference = step1_dependents.path + "/NC000962.fasta"
            print("Reference being used: %s" % reference)
            self.hqs = step1_dependents.path + "/15-3162-highqualitysnps.vcf"
            self.gbk_file = step1_dependents.path + "/NC_000962.gbk"
            self.upload_to = step1_dependents.upload_to
            self.remote = step1_dependents.remote
            self.script_dependents = step1_dependents.path
            
        if species == "para":
            self.dependents_dir="/mycobacterium/avium_complex/para_cattle-bison/script_dependents/script1"
            self.step1_dependents = Update_Directory(dependents_dir) #OBJECT not variable made
            self.spoligo_db = step1_dependents.path + "/nospoligo.txt"
            self.reference = step1_dependents.path + "/NC_002944.fasta"
            print("Reference being used: %s" % reference)
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
import os
import shutil
import sys

# vSNP classes
from vFunctions import Update_Directory
from vFunctions import Bruc_Private_Codes
from vFunctions import Get_Filters

class Step1:

    '''Instantiation, from vProperties import Step1, x = Step1("virus"), x.s'''

    def __init__(self, species):
        if species == "ab1":
            dependents_dir="/brucella/abortus1/script_dependents/script1"
            step1_dependents = Update_Directory(dependents_dir)

            spoligo_db = step1_dependents.path + "/nospoligo.txt"
            reference = step1_dependents.path + "/NC_00693c.fasta"
            print("Reference being used: %s" % reference)
            hqs = step1_dependents.path + "/NC_00693cHighestQualitySNPs.vcf"
            gbk_file = step1_dependents.path + "/NC_006932-NC_006933.gbk"
            
            self.option_list = [dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        
        if species == "ab3":
            dependents_dir="/brucella/abortus3/script_dependents/script1"
            step1_dependents = Update_Directory(dependents_dir)
            spoligo_db = step1_dependents.path + "/nospoligo.txt"
            reference = step1_dependents.path + "/CP007682-7683c.fasta"
            print("Reference being used: %s" % reference)
            hqs = step1_dependents.path + "/CP007682-7683cHighestQualitySNPs.vcf"
            gbk_file = step1_dependents.path + "/CP007682-CP007683.gbk"
                        
            self.option_list = [dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            
        if species == "canis":
            dependents_dir="/brucella/canis/script_dependents/script1"
            step1_dependents = Update_Directory(dependents_dir)
            spoligo_db = step1_dependents.path + "/nospoligo.txt"
            reference = step1_dependents.path + "/BcanisATCC23365.fasta"
            print("Reference being used: %s" % reference)
            hqs = step1_dependents.path + "/canisHighestQualitySNPs.vcf"
            gbk_file = step1_dependents.path + "/NC_010103-NC_010104.gbk"
                        
            self.option_list = [dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        
        if species == "ceti1":
            dependents_dir="/brucella/ceti1/script_dependents/script1"
            step1_dependents = Update_Directory(dependents_dir)
            spoligo_db = step1_dependents.path + "/nospoligo.txt"
            reference = step1_dependents.path + "/Bceti1Cudo.fasta"
            print("Reference being used: %s" % reference)
            hqs = step1_dependents.path + "/ceti1HighestQualitySNPs.vcf"
            gbk_file = None #script_dependents + "/no.gff"
                        
            self.option_list = [dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]

        if species == "ceti2":
            dependents_dir="/brucella/ceti2/script_dependents/script1"
            step1_dependents = Update_Directory(dependents_dir)
            spoligo_db = step1_dependents.path + "/nospoligo.txt"
            reference = step1_dependents.path + "/Bceti2-TE10759.fasta"
            print("Reference being used: %s" % reference)
            hqs = step1_dependents.path + "/ceti2HighestQualitySNPs.vcf"
            gbk_file = step1_dependents.path + "/NC_022905-NC_022906.gbk"
                        
            self.option_list = [dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]

        if species == "mel1":
            dependents_dir="/brucella/melitensis-bv1/script_dependents/script1"
            step1_dependents = Update_Directory(dependents_dir)
            spoligo_db = step1_dependents.path + "/nospoligo.txt"
            reference = step1_dependents.path + "/mel-bv1-NC003317.fasta"
            print("Reference being used: %s" % reference)
            hqs = step1_dependents.path + "/mel-bv1-NC003317-highqualitysnps.vcf"
            gbk_file = step1_dependents.path + "/NC_003317-NC_003318.gbk"
                        
            self.option_list = [dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        
        if species == "mel1b":
            dependents_dir="/brucella/melitensis-bv1b/script_dependents/script1"
            step1_dependents = Update_Directory(dependents_dir)
            spoligo_db = step1_dependents.path + "/nospoligo.txt"
            reference = step1_dependents.path + "/mel-bv1b-CP018508.fasta"
            print("Reference being used: %s" % reference)
            hqs = step1_dependents.path + "/mel-bv1b-CP018508-highqualitysnps.vcf"
            gbk_file = step1_dependents.path + "/mel-bv1b-CP018508.gbk"
                        
            self.option_list = [dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        
        if species == "mel2":
            dependents_dir="/brucella/melitensis-bv2/script_dependents/script1"
            step1_dependents = Update_Directory(dependents_dir)
            spoligo_db = step1_dependents.path + "/nospoligo.txt"
            reference = step1_dependents.path + "/mel-bv2-NC012441.fasta"
            print("Reference being used: %s" % reference)
            hqs = step1_dependents.path + "/mel-bv2-NC012441-highqualitysnps.vcf"
            gbk_file = step1_dependents.path + "/NC_012441-NC_012442.gbk"
                        
            self.option_list = [dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            
        if species == "mel3":
            dependents_dir="/brucella/melitensis-bv3/script_dependents/script1"
            step1_dependents = Update_Directory(dependents_dir)
            spoligo_db = step1_dependents.path + "/nospoligo.txt"
            reference = step1_dependents.path + "/mel-bv3-NZCP007760.fasta"
            print("Reference being used: %s" % reference)
            hqs = step1_dependents.path + "/mel-bv3-NZCP007760-highqualitysnps.vcf"
            gbk_file = step1_dependents.path + "/NZ_CP007760-NZ_CP007761.gbk"
                        
            self.option_list = [dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        
        if species == "suis1":
            dependents_dir="/brucella/suis1/script_dependents/script1"
            step1_dependents = Update_Directory(dependents_dir)
            spoligo_db = step1_dependents.path + "/nospoligo.txt"
            reference = step1_dependents.path + "/NC_017251-NC_017250.fasta"
            print("Reference being used: %s" % reference)
            hqs = step1_dependents.path + "/B00-0468-highqualitysnps.vcf"
            gbk_file = step1_dependents.path + "/NC_017251-NC_017250.gbk"
                        
            self.option_list = [dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            
        if species == "suis3":
            dependents_dir="/brucella/suis3/script_dependents/script1"
            step1_dependents = Update_Directory(dependents_dir)
            spoligo_db = step1_dependents.path + "/nospoligo.txt"
            reference = step1_dependents.path + "/NZ_CP007719-NZ_CP007718.fasta"
            print("Reference being used: %s" % reference)
            hqs = step1_dependents.path + "/highqualitysnps.vcf"
            gbk_file = None #script_dependents + "/no.gff"
                        
            self.option_list = [dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            
        if species == "suis4":
            dependents_dir="/brucella/suis4/script_dependents/script1"
            step1_dependents = Update_Directory(dependents_dir)
            spoligo_db = step1_dependents.path + "/nospoligo.txt"
            reference = step1_dependents.path + "/B-REF-BS4-40.fasta"
            print("Reference being used: %s" % reference)
            hqs = step1_dependents.path + "/suis4HighestQualitySNPs.vcf"
            gbk_file = None #script_dependents + "/no.gff"
       
        if species == "ovis":
            dependents_dir="/brucella/ovis/script_dependents/script1"
            step1_dependents = Update_Directory(dependents_dir)
            spoligo_db = step1_dependents.path + "/nospoligo.txt"
            reference = step1_dependents.path + "/BovisATCC25840.fasta"
            print("Reference being used: %s" % reference)
            hqs = step1_dependents.path + "/BovisATCC25840HighestQualitySNPs.vcf"
            gbk_file = step1_dependents.path + "/NC_009505-NC_009504.gbk"
            email_list = "tod.p.stuber@aphis.usda.gov"

            self.option_list = [dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            
        if species == "neo":
            dependents_dir="/brucella/neotomae/script_dependents/script1"
            step1_dependents = Update_Directory(dependents_dir)
            spoligo_db = step1_dependents.path + "/nospoligo.txt"
            reference = step1_dependents.path + "/KN046827.fasta"
            print("Reference being used: %s" % reference)
            hqs = step1_dependents.path + "/ERR1845155-highqualitysnps.vcf"
            gbk_file = step1_dependents.path + "/KN046827.gbk"
            email_list = "tod.p.stuber@aphis.usda.gov"

            self.option_list = [dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            
        if species == "af":
            dependents_dir="/mycobacterium/tbc/af2122/script_dependents/script1"
            step1_dependents = Update_Directory(dependents_dir)
            spoligo_db = step1_dependents.path + "/spoligotype_db.txt"
            reference = step1_dependents.path + "/NC_002945v4.fasta"
            print("Reference being used: %s" % reference)
            hqs = step1_dependents.path + "/highqualitysnps.vcf"
            gbk_file = step1_dependents.path + "/NC_002945v4.gbk"
                        
            self.option_list = [dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        
        if species == "h37":
            dependents_dir="/mycobacterium/tbc/h37/script_dependents/script1"
            step1_dependents = Update_Directory(dependents_dir)
            spoligo_db = step1_dependents.path + "/spoligotype_db.txt"
            reference = step1_dependents.path + "/NC000962.fasta"
            print("Reference being used: %s" % reference)
            hqs = step1_dependents.path + "/15-3162-highqualitysnps.vcf"
            gbk_file = step1_dependents.path + "/NC_000962.gbk"
                        
            self.option_list = [dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            
        if species == "para":
            dependents_dir="/mycobacterium/avium_complex/para_cattle-bison/script_dependents/script1"
            step1_dependents = Update_Directory(dependents_dir)
            spoligo_db = step1_dependents.path + "/nospoligo.txt"
            reference = step1_dependents.path + "/NC_002944.fasta"
            print("Reference being used: %s" % reference)
            hqs = step1_dependents.path + "/HQ-NC002944.vcf"
            gbk_file = step1_dependents.path + "/NC_002944.gbk"
                        
            self.option_list = [dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            
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
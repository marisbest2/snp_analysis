import os
import shutil

class Step1:

    '''Instantiation, from vProperties import Step1, x = Step1("virus"), x.s'''

    def __init__(self, species):

        def update_directory(dependents_dir): # UPDATE DIRECTORIES
            home = os.path.expanduser("~")
            print("dependents_dir %s\n" % dependents_dir)
            
            if os.path.isdir("/bioinfo11/TStuber/Results"): #check bioinfo from server
                upload_to = "/bioinfo11/TStuber/Results"
                remote="/bioinfo11/TStuber/Results" + dependents_dir
                if os.path.isdir("/Users/Shared"):
                    dep_path = "/Users/Shared"
                    dir_split = dependents_dir.split('/')[1:]
                    for i in dir_split:
                        dep_path += '/' + i
                        if not os.path.exists(dep_path):
                            os.makedirs(dep_path)
                    local = "/Users/Shared" + dependents_dir
                    if os.path.isdir(local):
                        try:
                            shutil.rmtree(local)
                            shutil.copytree(remote, local)
                        except:
                            pass
                elif os.path.isdir("/home/shared"):
                    dep_path = "/home/shared"
                    dir_split = dependents_dir.split('/')[1:]
                    for i in dir_split:
                        dep_path += '/' + i
                        if not os.path.exists(dep_path):
                            os.makedirs(dep_path)
                    local = "/home/shared" + dependents_dir
                    if os.path.isdir(local):
                        try:
                            shutil.rmtree(local)
                            shutil.copytree(remote, local)
                        except:
                            pass

        if species == "salmonella":
            dependents_dir="/gen-bact/salmonella/snp_pipeline/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/nospoligo.txt"
            reference = script_dependents + "/NC_016856-NC_016855.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/NC_016856-NC_016855HighestQualitySNPs.vcf"
            gbk_file = script_dependents + "/NC_016856-NC_016855.gbk"
            email_list = "tod.p.stuber@aphis.usda.gov"
            
            option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            return option_list, found
    
        if species == "ab1":
            dependents_dir="/brucella/abortus1/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/nospoligo.txt"
            reference = script_dependents + "/NC_00693c.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/NC_00693cHighestQualitySNPs.vcf"
            gbk_file = script_dependents + "/NC_006932-NC_006933.gbk"
            email_list = "tod.p.stuber@aphis.usda.gov"
            
            option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            return option_list, found
        
        if species == "ab3":
            dependents_dir="/brucella/abortus3/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/nospoligo.txt"
            reference = script_dependents + "/CP007682-7683c.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/CP007682-7683cHighestQualitySNPs.vcf"
            gbk_file = script_dependents + "/CP007682-CP007683.gbk"
            email_list = "tod.p.stuber@aphis.usda.gov"
            

            option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            return option_list, found
        
        if species == "canis":
            dependents_dir="/brucella/canis/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/nospoligo.txt"
            reference = script_dependents + "/BcanisATCC23365.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/canisHighestQualitySNPs.vcf"
            gbk_file = script_dependents + "/NC_010103-NC_010104.gbk"
            email_list = "tod.p.stuber@aphis.usda.gov"
            
            option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            return option_list, found

        if species == "ceti1":
            dependents_dir="/brucella/ceti1/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/nospoligo.txt"
            reference = script_dependents + "/Bceti1Cudo.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/ceti1HighestQualitySNPs.vcf"
            gbk_file = None #script_dependents + "/no.gff"
            email_list = "tod.p.stuber@aphis.usda.gov"
            
            option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            return option_list, found

        if species == "ceti2":
            dependents_dir="/brucella/ceti2/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/nospoligo.txt"
            reference = script_dependents + "/Bceti2-TE10759.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/ceti2HighestQualitySNPs.vcf"
            gbk_file = script_dependents + "/NC_022905-NC_022906.gbk"
            email_list = "tod.p.stuber@aphis.usda.gov"
            
            option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            return option_list, found

        if species == "mel1":
            dependents_dir="/brucella/melitensis-bv1/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/nospoligo.txt"
            reference = script_dependents + "/mel-bv1-NC003317.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/mel-bv1-NC003317-highqualitysnps.vcf"
            gbk_file = script_dependents + "/NC_003317-NC_003318.gbk"
            email_list = "tod.p.stuber@aphis.usda.gov"
            
            option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            return option_list, found

        if species == "mel1b":
            dependents_dir="/brucella/melitensis-bv1b/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/nospoligo.txt"
            reference = script_dependents + "/mel-bv1b-CP018508.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/mel-bv1b-CP018508-highqualitysnps.vcf"
            gbk_file = script_dependents + "/mel-bv1b-CP018508.gbk"
            email_list = "tod.p.stuber@aphis.usda.gov"
            
            option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            return option_list, found

        if species == "mel2":
            dependents_dir="/brucella/melitensis-bv2/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/nospoligo.txt"
            reference = script_dependents + "/mel-bv2-NC012441.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/mel-bv2-NC012441-highqualitysnps.vcf"
            gbk_file = script_dependents + "/NC_012441-NC_012442.gbk"
            email_list = "tod.p.stuber@aphis.usda.gov"
            
            option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            return option_list, found
            
        if species == "mel3":
            dependents_dir="/brucella/melitensis-bv3/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/nospoligo.txt"
            reference = script_dependents + "/mel-bv3-NZCP007760.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/mel-bv3-NZCP007760-highqualitysnps.vcf"
            gbk_file = script_dependents + "/NZ_CP007760-NZ_CP007761.gbk"
            email_list = "tod.p.stuber@aphis.usda.gov"
            
            option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            return option_list, found

        if species == "suis1":
            dependents_dir="/brucella/suis1/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/nospoligo.txt"
            reference = script_dependents + "/NC_017251-NC_017250.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/B00-0468-highqualitysnps.vcf"
            gbk_file = script_dependents + "/NC_017251-NC_017250.gbk"
            email_list = "tod.p.stuber@aphis.usda.gov"
            
            option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            return option_list, found

        if species == "suis3":
            dependents_dir="/brucella/suis3/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/nospoligo.txt"
            reference = script_dependents + "/NZ_CP007719-NZ_CP007718.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/highqualitysnps.vcf"
            gbk_file = None #script_dependents + "/no.gff"
            email_list = "tod.p.stuber@aphis.usda.gov"
            
            option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            return option_list, found

        if species == "suis4":
            dependents_dir="/brucella/suis4/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/nospoligo.txt"
            reference = script_dependents + "/B-REF-BS4-40.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/suis4HighestQualitySNPs.vcf"
            gbk_file = None #script_dependents + "/no.gff"

            email_list = "tod.p.stuber@aphis.usda.gov"
            
        if species == "ovis":
            dependents_dir="/brucella/ovis/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/nospoligo.txt"
            reference = script_dependents + "/BovisATCC25840.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/BovisATCC25840HighestQualitySNPs.vcf"
            gbk_file = script_dependents + "/NC_009505-NC_009504.gbk"
            email_list = "tod.p.stuber@aphis.usda.gov"

            option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            return option_list, found
            
        if species == "neo":
            dependents_dir="/brucella/neotomae/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/nospoligo.txt"
            reference = script_dependents + "/KN046827.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/ERR1845155-highqualitysnps.vcf"
            gbk_file = script_dependents + "/KN046827.gbk"
            email_list = "tod.p.stuber@aphis.usda.gov"

            option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            return option_list, found
            
        if species == "af":
            dependents_dir="/mycobacterium/tbc/af2122/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/spoligotype_db.txt"
            reference = script_dependents + "/NC_002945v4.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/highqualitysnps.vcf"
            gbk_file = script_dependents + "/NC_002945v4.gbk"
            email_list = "tod.p.stuber@aphis.usda.gov"
            
            option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            return option_list, found

        if species == "h37":
            dependents_dir="/mycobacterium/tbc/h37/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/spoligotype_db.txt"
            reference = script_dependents + "/NC000962.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/15-3162-highqualitysnps.vcf"
            gbk_file = script_dependents + "/NC_000962.gbk"
            email_list = "tod.p.stuber@aphis.usda.gov"
            
            option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            return option_list, found

        if species == "para":
            dependents_dir="/mycobacterium/avium_complex/para_cattle-bison/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/nospoligo.txt"
            reference = script_dependents + "/NC_002944.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/HQ-NC002944.vcf"
            gbk_file = script_dependents + "/NC_002944.gbk"
            email_list = "tod.p.stuber@aphis.usda.gov"
            
            option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            return option_list, found
class Step2:
    def __init__(self, species):
        if species == "virus":
            self.s = "virus"
        elif species == "bacteria":
            self.s= "bacteria"
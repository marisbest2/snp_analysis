import os
import sys
import subprocess
import glob
import re
import regex
import shutil
import time
import csv
import xlrd
import vcf
import xlsxwriter
import gzip
import zipfile
from numpy import mean
from collections import Counter
from datetime import datetime
from optparse import OptionParser
from concurrent import futures
from functools import partial
from collections import OrderedDict
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO
import git

###############################################


###############################################

class script1:

    '''A script1 object.  These objects have the following properties:
    Attributes:
        fastq_R1
        fastq_R2 '''

    def __init__(self, fastq_R1, fastq_R2, species=None):
        '''instantiate object whose name *name*.  Proved species to force refernce type'''
#        self.fastq_R1=fastq_R1
#        self.fastq_R2=fastq_R2

        global R1
        global R2
        R1 = fastq_R1
        R2 = fastq_R2
        
        global gbk_file
        gbk_file = ""
        global in_annotation_as_dict
        #in_annotation_as_dict = {}
        
        global directory
        global zips
        directory = str(os.getcwd())
        zips = directory + "/zips"
        
        os.makedirs(zips)
        shutil.move(R1, zips)
        shutil.move(R2, zips)
        
        R1 = zips + "/" + R1
        R2 = zips + "/" + R2
        
        global R1unzip
        global R2unzip
        R1unzip = re.sub('\.gz', '', R1)
        R2unzip = re.sub('\.gz', '', R2)

        self.species = species
        
        if self.species:
            species_force = species
            print("Sample will be ran as %s" % self.species)
        else:
            species_force=None

        global dependents_dir
        global reference
        global hqs
        global email_list
        global upload_to
        global remote
        global script_dependents
        global spoligo_db
        
        global mlst_type
        global octalcode
        global sbcode
        global binarycode
        global ave_coverage
        global zero_coverage_vcf
        global genome_coverage

        ### SET PARAMETERS
        if species_force:
            option_list, found = script1.parameters(species_force)
            dependents_dir = option_list[0]
            reference = option_list[1]
            hqs = option_list[2]
            gbk_file = option_list[3]
            email_list = option_list[4]
            upload_to = option_list[5]
            remote = option_list[6]
            script_dependents = option_list[7]
            spoligo_db = option_list[8]

        else:
            best_ref_found = script1.best_reference(self)
            self.species = best_ref_found
            try: # exit if best ref isn't in parameter list
                option_list, found = script1.parameters(best_ref_found)
                dependents_dir = option_list[0]
                reference = option_list[1]
                hqs = option_list[2]
                gbk_file = option_list[3]
                email_list = option_list[4]
                upload_to = option_list[5]
                remote = option_list[6]
                script_dependents = option_list[7]
                spoligo_db = option_list[8]
            except UnboundLocalError:
                print("\nprinted options list:")
                for i in option_list:
                    print(i)
                print("\nNo parameters options for best_ref_found or species given - UnboundLocalError\n\n")
                self.species = "NO FINDINGS"
            except TypeError:
                print("No parameters options for best_ref_found or species given - TypeError")
                self.species = "NO FINDINGS"
        try:
            shutil.copy2(reference, directory)
            shutil.copy2(hqs, directory)
        except FileNotFoundError:
            time.sleep(5)
            shutil.copy2(reference, directory)
            shutil.copy2(hqs, directory)
        except NameError:
            print("No parameters options for best_ref_found or species given - NameError")
            self.species = "NO FINDINGS"
        ### SET PARAMETERS

###############################################

    def show_fastqs(self):
        '''show FASTQs being used'''
        print("R1: %s\nR2: %s" % (R1, R2))
    
###############################################

    def sizeof_fmt(num, suffix='B'):
        for unit in ['','K','M','G','T','P','E','Z']:
            if abs(num) < 1024.0:
                return "%3.1f%s%s" % (num, unit, suffix)
            num /= 1024.0
        return "%.1f%s%s" % (num, 'Yi', suffix)
    
###############################################

    def unzipfiles(): # don't "self" if called from within
        for zip_filename in R1, R2:
            print("Unzipping... %s" % zip_filename)
            name_nogz = os.path.splitext(zip_filename)[0]
            write_out = open(name_nogz, 'wt')
            with gzip.open(zip_filename, 'rt') as f:
                file_content = f.read()
                print(file_content, file=write_out)
            write_out.close()

###############################################

    def update_directory(dependents_dir): # UPDATE DIRECTORIES
        home = os.path.expanduser("~")
        
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
                        try:
                            os.makedirs(dep_path)
                        except:
                            pass
                local = "/home/shared" + dependents_dir
                if os.path.isdir(local):
                    try:
                        shutil.rmtree(local)
                        shutil.copytree(remote, local)
                    except:
                        pass

        elif os.path.isdir("/Volumes/root/TStuber/Results"): #check bioinfo from Mac
            upload_to = "/Volumes/root/TStuber/Results"
            remote="/Volumes/root/TStuber/Results" + dependents_dir
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
                        try:
                            os.makedirs(dep_path)
                        except:
                            pass
                local = "/home/shared" + dependents_dir
                if os.path.isdir(local):
                    try:
                        shutil.rmtree(local)
                        shutil.copytree(remote, local)
                    except:
                        pass

        #### PLACE A CHECK FROM GITHUB
        
        elif os.path.isdir("/Users/Shared" + dependents_dir): #check local copy in shared folder
            upload_to ="not_found"
            remote = "not_found"
            local = "/Users/Shared" + dependents_dir
                
        elif os.path.isdir("/home/shared" + dependents_dir): #check local copy in shared folder
            upload_to ="not_found"
            remote = "not_found"
            local = "/home/shared" + dependents_dir
        
        elif os.path.isdir(home + "/dependencies" + dependents_dir): #check local copy from Git repo
            upload_to ="not_found"
            remote = "no remote"
            script_location = home # points to home directory
            local = home + "/dependencies" + dependents_dir # sets dependencies directory to home directory
        else:
            try:
                os.makedirs(home + "/dependencies")
                print("\n\nDOWNLOADING DEPENDENCIES FROM GITHUB... ***\n\n")
                git.Repo.clone_from("https://github.com/stuber/dependencies.git", home + "/dependencies")
                upload_to ="not_found"
                remote = "no remote"
                script_location = home # points to home directory
                local = home + "/dependencies" + dependents_dir # sets dependencies directory to home directory
            except FileExistsError:
                upload_to ="not_found"
                remote = "no remote"
                script_location = home # points to home directory
                local = home + "/dependencies" + dependents_dir # sets dependencies directory to home directory
                pass
            
        print("\n####################DIRECTORY LOCATION")
        print("####################upload_to: %s" % upload_to)
        print("####################remote: %s" % remote)
        print("####################local: %s\n" % local)
        
        return upload_to, remote, local

###############################################

    def get_annotations(line, in_annotation_as_dict): #for line in vfile
        #pos_found = False
        line=line.rstrip()
        if line.startswith("#"): # save headers to file
            return(line)
        elif not line.startswith("#"): # position rows
            #pos_found = False
            split_line = line.split('\t')
            chrom = split_line[0]
            position = split_line[1] # get position
        #print ("Getting annotations")
            for each_key, each_value in in_annotation_as_dict.items():
                pos_found = False
                if chrom == each_key:
                    for feature in each_value.features:
                        position = int(position)
                        #print(position)
                        if position in feature and "CDS" in feature.type:
                            myproduct = "none list"
                            mylocus = "none list"
                            mygene = "none list"
                            for p in feature.qualifiers['product']:
                                myproduct = p
                            for l in feature.qualifiers['locus_tag']:
                                mylocus = l
                            if "gene" in feature.qualifiers:
                                gene = feature.qualifiers['gene']
                                for g in gene:
                                    mygene = g
                            annotation_found = myproduct + ", gene: " + mygene + ", locus_tag: " + mylocus
                            pos_found = True
                if pos_found == False:
                    annotation_found = "No annotated product"
                #print(annotation_found)
                split_line[2] = annotation_found
                annotated_line = "\t".join(split_line)
                return(annotated_line)

###############################################

    def parameters(give_option):
        if give_option == "ab1":
            found=True
            #Remove network path at and left of "Results"
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
        
        if give_option == "ab3":
            found=True
            #Remove network path at and left of "Results"
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
        
        if give_option == "canis":
            found=True
            #Remove network path at and left of "Results"
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

        if give_option == "ceti1":
            found=True
            #Remove network path at and left of "Results"
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

        if give_option == "ceti2":
            found=True
            #Remove network path at and left of "Results"
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

        if give_option == "mel1":
            found=True
            #Remove network path at and left of "Results"
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

        if give_option == "mel2":
            found=True
            #Remove network path at and left of "Results"
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
            
        if give_option == "mel3":
            found=True
            #Remove network path at and left of "Results"
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

        if give_option == "suis1":
            found=True
            #Remove network path at and left of "Results"
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

        if give_option == "suis4":
            found=True
            #Remove network path at and left of "Results"
            dependents_dir="/brucella/suis4/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/nospoligo.txt"
            reference = script_dependents + "/B-REF-BS4-40.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/suis4HighestQualitySNPs.vcf"
            gbk_file = None #script_dependents + "/no.gff"

            email_list = "tod.p.stuber@aphis.usda.gov"
            
        if give_option == "ovis":
            found=True
            #Remove network path at and left of "Results"
            dependents_dir="/brucella/ovis/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/nospoligo.txt"
            reference = script_dependents + "/BovisATCC25840.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/BovisATCC25840HighestQualitySNPs.vcf"
            gbk_file = None #script_dependents + "/no.gff"
            email_list = "tod.p.stuber@aphis.usda.gov"

            option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            return option_list, found

        if give_option == "bovis":
            found=True
            #Remove network path at and left of "Results"
            dependents_dir="/mycobacterium/tbc/tbbov/script_dependents/script1"
            upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
            
            spoligo_db = script_dependents + "/spoligotype_db.txt"
            reference = script_dependents + "/NC_002945.fasta"
            print("Reference being used: %s" % reference)
            hqs = script_dependents + "/HighestQualitySNPs.vcf"
            gbk_file = script_dependents + "/NC_002945.gbk"
            email_list = "tod.p.stuber@aphis.usda.gov"
            
            option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
            return option_list, found
            
        if give_option == "af":
            found=True
            #Remove network path at and left of "Results"
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

        if give_option == "h37":
            found=True
            #Remove network path at and left of "Results"
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

        if give_option == "para":
            found=True
            #Remove network path at and left of "Results"
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
###############################################

    def finding_best_ref(v):
        count=0
        for fastq in R1unzip, R2unzip:
            with open(fastq) as in_handle:
                # all 3, title and seq and qual, were needed
                for title, seq, qual in FastqGeneralIterator(in_handle):
                    count += seq.count(v)
        return(v, count)

# FASTER TO UNZIP AND DO ABOVE THEN TO DO BELOW

#    def finding_best_ref(v):
#        count=0
#        for i in R1, R2:
#            with gzip.open(i, "rt") as handle:
#                for title, seq, qual in FastqGeneralIterator(handle):
#                    count += seq.count(v)
#        return(v, count)
###############################################
        
    def mlst(self):
    
        fastqs = glob.glob(zips + '/*.fastq')
        if len(fastqs) < 2:
            script1.unzipfiles()
        fastqs = glob.glob(zips + '/*.fastq')
        
        #https://bmcmicrobiol.biomedcentral.com/articles/10.1186/1471-2180-7-34
        write_ref = open("ST1-MLST.fasta", 'w')
        print(">ST1-MLST", file=write_ref)
        print("CGTTTCCCCAAGGAAGTGGAGGTTGCAGGCGATACGATCGATGTTGGCTACGGCCCGATCAAGGTTCATGCCGTCCGCAACCCGGCCGAACTGCCGTGGAAGGAAGAAAACGTCGATATCGCCCTTGAATGCACCGGCATTTTCACCTCGCGCGACAAGGCAGCACTTCATCTTGAAGCTGGCGCCAAGCGCGTCATCGTCTCCGCTCCCGCAGACGGTGCCGATCTCACCGTCGTCTATGGTGTCAACAACGACAAGCTGACGAAGGACCATCTGGTCATCTCCAACGCTTCGTGTACCACCAACTGCCTTGCGCCGGTGGCTCAGGTTCTCAACGATACTATCGGTATCGAAAAGGGCTTCATGACCACGATCCACTCCTATACGGGCGACCAGCCGACGCTGGACACCATGCACAAGGATCTCTACCGCGCCCGCGCCGCTGCCCTTTCCATGATCCCGACCTCGACGGGTGCGGCCAAGGCCGTCGGTCTCGTTCTGCCGGAACTGAAAGGCAAGCTCGACGGCGTTGCCATTCGCGTCCCGACCCCAAATGTCTCGGTCGTTGATCTCACCTTCATCGCCAAGCGTGAAACCACCGTTGAAGAAGTCAACAATGCGATCCGCGAAGCCGCCAATGGCCGCCTCAAGGGCATTCTCGGCTATACCGATGAGAAGCTCGTCTCGCACGACTTCAACCACGATTCCCATTCCTCGGTCTTCCACACCGACCAGACCAAGGTTATGGACGGCACCATGGTGCGTATCCTGTCGTGGTACGACAATGAATGGGGCTTCTCCAGCCGCATGAGCGACACCGCCGTCGCTTTGGGCAAGCTGATCTGATAACGGCAACGCCTCTCCTTCACTGGCGAGGCGTTTTCATTTCTTGATAAGGACCGAGAGAAGAAACATGATGTTCCGCACCCTTGACGATGCCAATGTCCAATCCAAGCGCGTGCTGGTCCGTGTTGACCTCAACGTGCCGAAAATCCGCCGTTCTGCTCGCTGGTCTTAACACCCCGGGCGTCACCACCGTGATCGAGCCGGTCATGACGCGCGATCATACGGAAAAGATGCTGCAAGACTTTGGCGCAGACCTGACGGTTGAAACCGATAAGGATGGTGTGCGCCATATCCGTATTGTCGGCCAGGGCAAGCTTACCGGCCAGACCATCGACGTGCCGGGTGATCCCTCGTCAACGGCTTTTCCGCTGGTGGCCGCCCTTCTGGTCGAAGGTTCGGAGGTCACCATCCGCAATGTGCTGATGAACCCGACCCGCACCGGCCTGATCCTGACGTTGCAGGAAATGGGGGCGGATATCGAGATCATCGATCCACGCCTTGCCGGCGGCGAGGATGTCGCCGATCTGCGCGTCAAGGCCTCGAAGCTGAAAGGCGTTGTCGTTCCGCCGGAACGTGCGCCTTCGATGATCGATGAATATCCGGTTCTGGCCATTGCCGCGTCTTTTGCGGAAGGCGAAACCGTGATGGACGGTCTCGATGAACTGCGCGTCAAGGAATCGGATCGTCTGGCGGCCGTTGCGCGCGGCCTTGAAGCCAATGGTGTCGATTGTACCGAAGGCGAGATGTCGCTGACGGTTCGTGGCCGCCCCGGCGGCAAGGGGCTGGGCGGTGGCACGGTTGCAACCCACCTCGACCACCGCATCGCGATGAGTTTCCTCGTCATGGGCCTTGCATCGGAAAAGCCGGTTACGGTGGATGACAGCACCATGATCGCCACCTCTTTCCCGGAATTCATGGGCATGATGGCGGGGCTGGGGGCGAAGATTGCCGAAAGCGGTGCAGAATGAAATCGTTCGTCGTCGCCCCGTTCATTGTCGCCATTGACGGACCGGCCGCCTCGGGCAAGGGAACCCTTGCCCGGCGGATCGCGACACATTACGGGATGCCGCATCTCGATACGGGCCTGACCTATCGCGCGGTCGCCAAGAGCCGCGCTCTGTCATTCTGGCCGTGGCAGGCCCGGTGGACGGCGACGAGATCGACCTCACCAATTGCGACTGGGTCGTGCGTCCTAAAAAGATGATCGCTGATCTGGGCTTTGAAGACGTGACCGTCCTCAATGATTTCGAGGCGCAGGCCCTTGCCGTGGTTTCGCTGGAAGGCCACCATATGGAACAGATCGGCGGCAAACCGGAGGAGGCTGTTGCCACCCGCGTCGTGCTCGGCCCCGGCACGGGCCTTGGCGTGGCAGGTCTGTTTCGCACACGTCATGCATGGGTTCCGGTTCCCGGTGAAGGCGGTCATATCGATATCGGTCCACGCACCGAACGCGACTACCAGATTTTCCCGCATATCGAACGCATCGAAGGGCGTGTCACCGGCGAGCAAATTCTTAGCGGGCGGGGCCTGCGCAACCTCTATCTGGGCATCTGCGCCGCCGACAAGATCACGCCCACCCTTGAGACGCCAGTAGACATTACATCCGCCGGACTGGACGGCAGCAATCCACAAGCCGCAGAAACGCTTGACCTCTTCGCCACCTATCTGGGGCGGCTTGCGGGCGACCTTGCGCTCATTTTCATGGCGCATGGCGGCGTTTATCTTTCGGGTGGCATCCCGGTGCGCATCCTTTCCGCCCTCAAGGCCGGTTCGTTCCGCGCAACCTTCGAGGACAAGGCCCCGCACAAGGCCATCATGCGCGACATACCGGTCCGCGTTATCACATATCAACTGGCGGCCTTAACCGGGCTTTCCGCTTTCGCCCGCACCCCCTCGCGCTTTGAAGTTTCGACCGAGGGCCGCCGCTGGCGCATGCGCCGCTAGAGCATTTCCGAGCCAAAAGTGCGAAGCGGTTCCGTTTCCCAACGAGCCGACCGCGGCTGCGCTTGCCTATGGTCTCGACAAGAGCGAAGGCAAGACCATCGCTGTCTATGACCTTGGCGGCGGTACTTTCGACGTGTCGGTTCTGGAAATCGGCGACGGCGTTTTTGAAGTGAAGTCCACCAATGGCGACACGTTCCTTGGCGGTGAAGACTTCGATATTCGTCTGGTCGAATATCTGGTTGCCGAGTTCAAGAAGGAAAGTGGCATCGACCTGAAGAACGACAAGCTTGCCCTGCAGCGCCTCAAGGAAGCTGCCGAAAAGGCCAAGATCGAACTGTCGTCCTCGCAGCAGACCGAAATCAACCTGCCGTTCATCACGGCTGACCAGACTGGCCCGAAGCATCTGGCGATCAAGCTGTCGCGCGCCAAGTTTGAAAGCCTGGTCGATGATCTCGTGCAGCGCACGGTCGAGCCGTGCAAGGCGGCGCTCAAGGATGCCGGCCTCAAGGCTGGCGAAATTGACGAAGTGGTTCTGGTCGGCGGCATGACCCGCATGCCCAAGATTCAGGAAGTCGTGAAGGCCTTCTTCGGCAAGGAACCGCACAAGGGCGTGAACCCGGATGAAGTCGTGGCCATGGGCGCGGCGATCCAGGGCGGCGTTTTGCAGGGCGACGTCAAGGACGTGCTGCTGCTCGACGTGACCCCGCTTTCGCTCGGCATTGAAACGCTGGGCGGCGTGTTCACCCGCCTGATCGAACGCAACACCACTATCCCGACCAAGAAGTCGCAGACCTTCTCCACGGCTGAGGACAACCAGTCGGCCGTGACGATCCGCGTCTTCCAGGGCGAGCGTGAAATGGCAGCCGATAACAAGCTGCTTGGACAGTTCGACCTCGTTGGCATTCCGCCACGTCCCTGCCCGGAAAGCTTGCCGATTGCCAGGAGCGCGATCCGGCCAAGTCCGAAATCTTCATCGTCGAGGGCGATTCGGCAGGCGGTTCCGCCAAGAGCGGGCGCTCGCGCCAGAATCAGGCCATTCTGCCGCTGCGCGGCAAAATCCTGAACGTGGAACGCGTGCGTTTCGACCGGATGATTTCATCCGATCAGGTGGGCACCCTCATCACGGCGCTTGGCACCTCCATCGGCAAGGATGAAACGCACGGCTTCAACGCCGACAAGCTGCGTTATCACAAGATCATCATCATGACCGACGCCGACGTCGATGGCGCCCATATTCGTACGCTTCTGCTCACCTTCTTCTTCCGGCAGATGCCGGAACTGATCGAACGCGGGCATATCTATATCGCGCAGCCGCCGCTCTATAAGGTGACACGCGGCAAGTCTTCGCAATATATCAAGAACGAAGCCGCCTTTGAGGATTTCCTCATCGAAACCGGCCTTGAAGAAACGACACTGGAACTGGTGACTGGCGAAATGCGCGCCGGGCCGGATTTGCGCTCGGTGGTGGAGGATGCGCGCACGCTGCGTCAGCTTCTGCACGGCCTGCACACCCGCTATGACCGCAGCGTGGTGGAACAGGCGGCAATTGCCGGCCTGCTCAACCCCGATGCCTCAAGGGACAATGCAACGGCACAGCATTCCGCCGATACGGTTGCCAAGCGTCTCGACATGATTTCGGAAGAGACCGAGCGCGGCTGGAGCGGCCATGTGATGGAAGACGGCGGCTATCGCTTCGAGCGTATGGTGCGCGGTGTAAAGGATATCGCCATTCTCGACATGGCCCTGCTCGGCTCGGCCGATGCCCGCCAGGTCGACCGAGATCGAGATGTATTCCCGCCTGATCCATACGGTCGATCATATCGAAGGCCGCCTGCGTGACGGCATGGATGCGTTTGACGGCTTCCTCAGCCATGCATGGGCTGTGACGGTGACAGGCGCGCCGAAGCTGTGGGCAATGCGCTTTCTTGAGGAAAACGAACGCAGCCCGCGCGCATGGTATGGCGGCGCGATCGGCATGATGCATTTCAATGGCGATATGAATACAGGGCTGACGCTGCGCACCATCCGCATCAAGGATGGTGTGGCGGAAATCCGTGCAGGGGCGACGCTTCTGTTCGATTCCAACCCTGACGAGGAAGAAGCCGAGACCGAATTGAAGGCATCGGCCATGATTGCGGCTGTGCGGGACGCACAGAAGAGCAATCAGATCGCGGAAGAAAGTGTGGCGGCAAAGGTGGGTGAGGGGGTTTCGATCCTGCTGGTCGATCACGAGGATTCCTTCGTCCATACGCTTGCCAATTATTTCCGCCAGACGGGCGCCAAGGTTTCCACCGTGCGTTCACCGGTGGCAGAGGAGATATTCGACCGCGTCAATCCCGATCTGGTGGTGTTATCGCCGGGACCGGGCTCGCCGCAGGATTTCGATTGCAAGGCGACCATCGATAAGGCGCGCAAGCGCCAGCTTCCGATTTTTGGCGTCTGCCTCGGCCTTCAGGCACTGGCGGAAGCCTATGGCGGGGCGTTGCGCCAGCTTCGCGTTCCGGTGCATGGCAAGCCTTCACGCATCCGCGTATCAAAGCCGGAGCGCATTTTCTCCGGCTTGCCGGAGGAAGTGACGGTGGGGCGTTATCATTCGATCTTCGCCGATCCTGAACGCCTGCCGGATGATTTTCTCGTCACAGCCGAAACGGAAGACGGGATCATAGCCTGCGGTGGAGGTGGTGATGGTGCCGCCGGGCTCCAGCCTGCCTGCGGATGCGGGGCTTGTCGTGTTGCCCGGCACCAAATCCACGATTGCCGATCTGCTGGCGCTGCGTGAAAACGGCTGGGACCGCGAATTGGTCGCCCATGTGAAGCGGGGCGGGCATGTGCTTGGTATTTGCGGCGGGTTTCAAATGCTTGGACGGCGGATCAGTGACCCGGCGGGTATTGAAGGCAATGTGCGCGATATCGAGGGGCTGGGCCTTCTCGATATCGAGACGATGACGGAGCCGGAAAAAGTGGTTCGCAATGTTGAGGCGGTGTCGCTGCTGCATGATGAGCCGCTGGAGGGCTATGAAATCCACATCGGGCGCACCAGCGGGCCGGATATGGCGCGGCCATTTGCGCGTATCGGCGATCATGATGATGGGGCCGTCTCGCCCGATGGTCGTATCATGGGAACCTATCTCCACGGTATTTTCAGTGCGGATCGTTTCCGCCACCACTTTTTGCGCGCGCTGGGTGTGGAAGGCGGCCAGATGAATTATCGCGAGAGCGTCGAAGAGGCTCTGGGCGAACTGGCTGAAGGGCTGGAAGCCTCGCTGGATATTGATGGCCTGTTTGCGCTGGCATGATTGACGCCGCGAAGCCGAAAGCCTAGTGTCAAACCATGTGACAGGTTTTGCCGGAACGAATCCCCGGCAATACCAAAAGGGAATGCGACGGACGGACCCACGCCGGGCGTCTTTATCGCAGCCGACCCCGCGACTGTAGAGCGGAGAGGGAAGAGGCAAGCCGGGCAACCGGCAGCCACTGGAAATCAGATGCGATAATGCAACATCGCATTTTTGCCATCTTCTCGACAGATTATCTCCACACAATGGGGCATTTCGTGCCGCAATTACCCTCGATATGTCACCCCTGTCAGCGCGGCATGGGCGGTTTACTCCCGATGCTGCCCGCCCGATAAGGGACCGCGCAAAACGTAATTTGTGTAAGGAGAATGCCATGCGCACTCTTAAGTCTCTCGTAATCGTCTCGGCTGCGCTGCTGCCGTTCTCTGCGACCGCTTTTGCTGCCGACGCCATCCAGGAACAGCCTCCGGTTCCGGCTCCGGTTGAAGTAGCTCCCCAGTATAGCTGGGCTGGTGGCTATACCGGTCTTTACCTTGGCTATGGCTGGAACAAGGCCAAGACCAGCACCGTTGGCAGCATCAAGCCTGACGATTGGAAGGCTGGCGCCTTTGCTGGCTGGAACTTCCAGCAGGACCAGATCGTATACGGTGTTGAAGGTGATGCAGGTTATTCCTGGGCCAAGAAGTCCAAGGACGGCCTGGAAGTCAAGCAGGGCTTTGAAGGCTCGCTGCGTGCCCGCGTCGGCTACGACCTGAACCCGGTTATGCCGTACCTCACGGCTGGTATTGCCGGTTCGCAGATCAAGCTTAACAACGGCTTGGACGACGAAAGCAAGTTCCGCGTGGGTTGGACGGCTGGTGCCGGTCTCGAAGCCAAGCTGACGGACAACATCCTCGGCCGCGTTGAGTACCGTTACACCCAGTACGGCAACAAGAACTATGATCTGGCCGGTACGACTGTTCGCAACAAGCTGGACACGCAGGATATCCGCGTCGGCATCGGCTACAAGTTCTAATTATAGCATAATTGGACACGGAAAACCGGACAGGCAACTGTCCGGTTTTTTGTTGTCTGCAAAGGTCGAGAAAGCGCGGCAGAGCAACGGCGGCAGCCTGATTTTCAGGGGAAATGAAGTGGAGGCTTCTGTTGCCAGGTGCCTCCGAACCCCGCCTTAAGGGGCTAACCCTAAGGACTTTAGAGTGGGTTTCCCGCACCGCCATTAGGCAGCGAGAGCATAACCCTGAGCATTGTTGTCATTTGCAACTACTCTGTTGACCCGATAACGGTGGTATCATGCCGAGTAAAAGAGCGATCTTTACACCCTTGTCGATCCTGTTTCGCCCCCGCCACAACACAGCCTGATCGGCAAGCTGTGCTGTGGTGGAGGCGCCGGGTACCGCCCCCGGGTCCAATGGGTTTATTACACCGTCCGTTTATCACCATAGTCGGCTTGCGCCGACAGGACGTATATAGGCGTGGTTTTTACCGATTGGAAGGGGGCTTGTGCGTTTTCGCGCAAGACCGACAGAGGTGGTGCGGCCCTTCCGTTCATTTTCCATTGACAGCTTCCGCGTGCTGGTCAATCCTCACAATATATCGGGATCGGCCTTGAAGAGGCTTGGCGCAGCCGGGGCGGAAACCATGGCTGAAACGGGGACGATATGCCCCAATCGAAGGAGAGTGGATATATGAGTGAATATCTCGCGGATGTCCGTCGCTATGATGCTGCCGCCGATGAGGCCGTTGTCGAGAAAATCGTCAAGCATCTTGGCATTGCGCTTCGCAATCGCGATTCCTCGCTCGTTTCGGCAAGC", file=write_ref)
        write_ref.close()

        directory = str(os.getcwd())
        print(directory)
        sample_reference_mlst_location = directory + "/ST1-MLST.fasta"
        read_base = os.path.basename(R1)
        sample_name = re.sub('_.*', '', read_base)
        sample_name_mlst = sample_name + "-mlst"
        print ("mlst reference: %s" % sample_reference_mlst_location)
        ref=re.sub('\.fasta', '', os.path.basename(sample_reference_mlst_location))
        print(ref)

        loc_sam_mlst=directory + "/" + sample_name_mlst
        print ("\n--")
        print("sample name mlst: %s" % sample_name_mlst)
        print("sample reference mlst: %s" % sample_reference_mlst_location)
        print("ref, no ext: %s " % ref)
        print("Read 1: %s" % R1)
        print("Read 2: %s\n" % R2)
        print("directory: %s" % directory)
        print("loc_sam_mlst: %s" % loc_sam_mlst)
        print ("--\n")

        os.system("samtools faidx {}" .format(sample_reference_mlst_location))
        os.system("picard CreateSequenceDictionary REFERENCE={} OUTPUT={}" .format(sample_reference_mlst_location, directory + "/" + ref + ".dict"))
        os.system("bwa index {}" .format(sample_reference_mlst_location))
        
        print("\n@@@ BWA mem")
        samfile_mlst = loc_sam_mlst + ".sam"
        os.system("bwa mem -M -t 16 {} {} {} > {}" .format(sample_reference_mlst_location, R1, R2, samfile_mlst))
        
        print("\nAdd read group and out all BAM")
        all_bam_mlst = loc_sam_mlst + "-all.bam"
        os.system("picard AddOrReplaceReadGroups INPUT={} OUTPUT={} RGLB=lib1 RGPU=unit1 RGSM={} RGPL=illumina" .format(samfile_mlst, all_bam_mlst, sample_name_mlst))

        print("\n@@@ Samtools mapped")
        mapbam = loc_sam_mlst + "-unmapped.bam"
        os.system("samtools view -h -F4 -b -T {} {} -o {}" .format(sample_reference_mlst_location, all_bam_mlst, mapbam))

        print("\n@@@ Sort BAM")
        sortedbam = loc_sam_mlst + "-sorted.bam"
        os.system("samtools sort {} -o {}" .format(mapbam, sortedbam))

        print("\n@@@ Index BAM")
        os.system("samtools index {}" .format(sortedbam))

        print("\n@@@ Calling SNPs with UnifiedGenotyper")
        vcf_mlst = directory + "/" + sample_name + "_mlst" + ".vcf"
        os.system("gatk -R {} -T UnifiedGenotyper -glm BOTH -out_mode EMIT_ALL_SITES -I {} -o {} -nct 8" .format(sample_reference_mlst_location, sortedbam, vcf_mlst))

        # Position 1629 was too close to the end of glk sequence.  Reads would not assemble properly to call possilbe SNP, therefore 100 bases of the gene were added.  Because of this all positions beyond this point are 100 more.  Same with position 1645 and 2693.

        target_vcf_positions = [231, 297, 363, 398, 429, 523, 631, 730, 1247, 1296, 1342, 1381, 1648, 1685, 1741, 1754, 2165, 2224, 2227, 2297, 2300, 2344, 2352, 2403, 2530, 2557, 2578, 2629, 3045, 3054, 3118, 3295, 3328, 3388, 3966, 3969, 4167, 4271, 4296, 4893, 4996, 4998, 5058, 5248, 5672, 5737, 5928, 5963, 5984, 5987, 6025, 6045, 6498, 6499, 6572, 6627, 6715, 6735, 6745, 6785, 6810, 6828, 6845, 6864, 6875, 7382, 7432, 7464, 7594, 7660, 7756]

        pos_call_dict={}
        vcf_reader = vcf.Reader(open(vcf_mlst, 'r'))
        for record in vcf_reader:
            if record.ALT[0]:
                pos_call_dict.update({record.POS:record.ALT[0]})
            else:
                pos_call_dict.update({record.POS:record.REF})

        mlst_string=[]
        for i in target_vcf_positions:
            if i in pos_call_dict:
                mlst_string.append(pos_call_dict[i]) #if number in list then get value
        mlst_join =  ''.join(map(str, mlst_string))
        print(mlst_join)

        mlst_dictionary={}
        mlst_dictionary["CTCCCGGGGCGACCCGATCGAAGCGGGAAGGCCACGGCGCGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA"] = "MLST type 01"
        mlst_dictionary["CTCCCGGGGCGACCCGAGCGAAGCGGGAAGGCCACGGCGCGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA"] = "MLST type 02"
        mlst_dictionary["CTCCCGTGGCGACCCGAGCGAAGCGGGAAGGCCACGGCGCGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA"] = "MLST type 03"
        mlst_dictionary["CTCCCGGGGCGACCCGAGCGAAGCGGGAAGGCCAAGGCGCGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA"] = "MLST type 04"
        mlst_dictionary["CTCCCGGGGCGACCCGATCGAAGCGGGAAGGCCACGGCGAGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA"] = "MLST type 05"
        mlst_dictionary["TTCCTGGGGCAACCCGAGCGAGGCAGGGAGGCCGCGGCTCGTGAGCGGTCGGGCATCTGTCCCGCGGGGTA"] = "MLST type 06"
        mlst_dictionary["CTTCCTGGCCGAGCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCTCGGCTGGGTACCTGTCTCGCGGTGCT"] = "MLST type 07"
        mlst_dictionary["CTTCCTGGCCGACCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCGCGGCTGGGTACCCGTCTCGCGGTGCT"] = "MLST type 08"
        mlst_dictionary["CTTCCTGGCCGACCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCGCGGCTGGGTACCTGTCTCGTGGTGCT"] = "MLST type 09"
        mlst_dictionary["CTTCCTGGCCGACCCGAGTGAAGGGGGGGGGCCACGGCGCGTGCTCGGCTGGGTACCTGTCTCGCGGTGCT"] = "MLST type 10"
        mlst_dictionary["CTTCCTGGCCGACCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCGCGGCTGGGTACCTGTCTCGCGGTGCT"] = "MLST type 11"
        mlst_dictionary["CTTCCTGGCCGACCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCTCGGCTGGGTACCTGTCTCGCGGTGCT"] = "MLST type 12"
        mlst_dictionary["CCCCCGGGCCGACTCGAGCGAAGCGAAGAGGCCACGGCGCGTGAGTGACCAGGCACCTATCCCACGGGGTA"] = "MLST type 13"
        mlst_dictionary["CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAGTGCGTGAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 14"
        mlst_dictionary["CCCCCGGGCCGACCCGGGCGAAGCGGGGAGGCTACGGTGCGTGAGTGGCCAGGCACCTGTCCCGCGAGGTA"] = "MLST type 15"
        mlst_dictionary["CCCCCGGCCCGACCCGGGCGAAGCGGGGAGGCTACGGTGCGTGAGTGGCCAGGCACCTGTCCCGCGAGGTA"] = "MLST type 16"
        mlst_dictionary["CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAATGCGTGAGTGGCCAGGCACCTGTCCCGCAGGGTA"] = "MLST type 17"
        mlst_dictionary["CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAATGCGTGAGTGGCCAGGCACCTGTCCCGCAGGCTA"] = "MLST type 18"
        mlst_dictionary["CCCCCGGGCCGACCCGAGCGAAGCGGGGAGGACACGGCGCGTGAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 19"
        mlst_dictionary["CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAATGCGTGAGTGGCCAGGCACATGTCCCGCAGGGTA"] = "MLST type 20"
        mlst_dictionary["CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAATGCGTGAGTGGCCAGGCACATGCCCCGCAGGGTA"] = "MLST type 21"
        mlst_dictionary["CCCCCGGGCCGACCCGAGCGAGGCGGGGAGGCCACGGCGCGGGAGTGGCCAGACACCTGTCCTGCGGGGTA"] = "MLST type 22"
        mlst_dictionary["CCCCCGGGCTGACCCGAGCGAAACGGGGAAGCCACGGCGCGTAAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 23"
        mlst_dictionary["CCCCCGGGCTGACCCGAGCGGAACGGGGAAGCCACGGCGCGTAAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 23x"
        mlst_dictionary["CCCCCGGGCCGACCCGAGCAAAGCGGGGAGGCCACGGCGCGTAAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 24"
        mlst_dictionary["CCCCCGGGCCGACCCGAGCGAAGCGGGGAGGCCACGGCGCGTAAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 25"
        mlst_dictionary["CCCCCGGGCCGACCCGAGCGAAGCGGGGAGGCCACGGCGCGTAAGTGGCCAAGCACCTGTTCCGCGGGGTA"] = "MLST type 26"
        mlst_dictionary["CCCCCGGGCCGACCCGAGCGAAGCGGGGAGACCACGGCGCATAAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 27"
        mlst_dictionary["CCCTCGGGCCGACCTGAGCGAAGCGGGGAGACCACGGCGCATAAGTGGCCAGGCTCCTGTCCCGCGGGGTA"] = "MLST type 28"

        for i in fastqs: #remove unzipped fastq files to save space
            os.remove(i)

        remove_files = glob.glob('ST1-MLST*')
        for i in remove_files:
            os.remove(i)
        remove_files = glob.glob('*-mlst*')
        for i in remove_files:
            os.remove(i)
        remove_files = glob.glob('*_mlst.vcf.idx')
        for i in remove_files:
            os.remove(i)

        write_out = open("mlst.txt", 'w')
        
        if mlst_join in mlst_dictionary:
            mlst_type = mlst_dictionary[mlst_join]
            print(mlst_type)
            print(mlst_type, file=write_out)
        else:
            print("NO MLST MATCH FOUND\n")
            print("NO MLST MATCH FOUND", file=write_out)
        
        write_out.close()

###############################################

    def finding_sp(v):
        total=0
        total_finds=0
        for fastq in R1unzip, R2unzip:
            with open(fastq) as in_handle:
                # all 3, title and seq and qual, were needed
                for title, seq, qual in FastqGeneralIterator(in_handle):
                    #if total < 6: # doesn't make a big different.  Might as well get full counts
#                    total += sum(seq.count(x) for x in (v)) #v=list of for and rev spacer
                    for spacer in v:
                        total_finds = len(regex.findall("(" + spacer + "){s<=1}", seq, overlapped=True))
                        total += total_finds
        return(v, total)
        #m = regex.findall("(bab){s<=1}", seq, overlapped=True)
###############################################

    def binary_to_octal(binary):
        binary_len = len(binary)
        i = 0
        ie = 1
        octal = ""
        while ie < 43:
            ie = i + 3
            print(binary[i:ie])
            region = binary[i:ie]
            region_len = len(region)
            i += 3
            if int(region[0]) == 1:
                if region_len < 2: # for the lone spacer 43.  When present needs to be 1 not 4.
                    oct = 1
                else:
                    oct = 4
            else:
                oct = 0
            try:
                if int(region[1]) == 1:
                    oct += 2
                if int(region[2]) == 1:
                    oct += 1
            except IndexError:
                pass
            octal = octal + str(oct)
        return(octal)
    
###############################################

    def binary_to_hex(binary):
        section1 = binary[0:7]
        section2 = binary[7:14]
        section3 = binary[14:21]
        section4 = binary[21:28]
        section5 = binary[28:36]
        section6 = binary[36:43]

        hex_section1 = hex(int(section1, 2))
        hex_section2 = hex(int(section2, 2))
        hex_section3 = hex(int(section3, 2))
        hex_section4 = hex(int(section4, 2))
        hex_section5 = hex(int(section5, 2))
        hex_section6 = hex(int(section6, 2))

        return(hex_section1.replace('0x', '').upper() + "-" + hex_section2.replace('0x', '').upper() + "-" + hex_section3.replace('0x', '').upper() + "-" + hex_section4.replace('0x', '').upper() + "-" + hex_section5.replace('0x', '').upper() + "-" + hex_section6.replace('0x', '').upper())

###############################################

    def spoligo(self):
        
        print("\nFinding spoligotype pattern...\n")
        
        fastqs = glob.glob(zips + '/*.fastq')
        if len(fastqs) < 2:
            script1.unzipfiles()
        fastqs = glob.glob(zips + '/*.fastq')
      
        '''spoligo spacers'''
        spoligo_dictionary = {}
        spoligo_dictionary["spacer01"] = ["TGATCCAGAGCCGGCGACCCTCTAT", "ATAGAGGGTCGCCGGCTCTGGATCA"]
        spoligo_dictionary["spacer02"] = ["CAAAAGCTGTCGCCCAAGCATGAGG", "CCTCATGCTTGGGCGACAGCTTTTG"]
        spoligo_dictionary["spacer03"] = ["CCGTGCTTCCAGTGATCGCCTTCTA", "TAGAAGGCGATCACTGGAAGCACGG"]
        spoligo_dictionary["spacer04"] = ["ACGTCATACGCCGACCAATCATCAG", "CTGATGATTGGTCGGCGTATGACGT"]
        spoligo_dictionary["spacer05"] = ["TTTTCTGACCACTTGTGCGGGATTA", "TAATCCCGCACAAGTGGTCAGAAAA"]
        spoligo_dictionary["spacer06"] = ["CGTCGTCATTTCCGGCTTCAATTTC", "GAAATTGAAGCCGGAAATGACGACG"]
        spoligo_dictionary["spacer07"] = ["GAGGAGAGCGAGTACTCGGGGCTGC", "GCAGCCCCGAGTACTCGCTCTCCTC"]
        spoligo_dictionary["spacer08"] = ["CGTGAAACCGCCCCCAGCCTCGCCG", "CGGCGAGGCTGGGGGCGGTTTCACG"]
        spoligo_dictionary["spacer09"] = ["ACTCGGAATCCCATGTGCTGACAGC", "GCTGTCAGCACATGGGATTCCGAGT"]
        spoligo_dictionary["spacer10"] = ["TCGACACCCGCTCTAGTTGACTTCC", "GGAAGTCAACTAGAGCGGGTGTCGA"]
        spoligo_dictionary["spacer11"] = ["GTGAGCAACGGCGGCGGCAACCTGG", "CCAGGTTGCCGCCGCCGTTGCTCAC"]
        spoligo_dictionary["spacer12"] = ["ATATCTGCTGCCCGCCCGGGGAGAT", "ATCTCCCCGGGCGGGCAGCAGATAT"]
        spoligo_dictionary["spacer13"] = ["GACCATCATTGCCATTCCCTCTCCC", "GGGAGAGGGAATGGCAATGATGGTC"]
        spoligo_dictionary["spacer14"] = ["GGTGTGATGCGGATGGTCGGCTCGG", "CCGAGCCGACCATCCGCATCACACC"]
        spoligo_dictionary["spacer15"] = ["CTTGAATAACGCGCAGTGAATTTCG", "CGAAATTCACTGCGCGTTATTCAAG"]
        spoligo_dictionary["spacer16"] = ["CGAGTTCCCGTCAGCGTCGTAAATC", "GATTTACGACGCTGACGGGAACTCG"]
        spoligo_dictionary["spacer17"] = ["GCGCCGGCCCGCGCGGATGACTCCG", "CGGAGTCATCCGCGCGGGCCGGCGC"]
        spoligo_dictionary["spacer18"] = ["CATGGACCCGGGCGAGCTGCAGATG", "CATCTGCAGCTCGCCCGGGTCCATG"]
        spoligo_dictionary["spacer19"] = ["TAACTGGCTTGGCGCTGATCCTGGT", "ACCAGGATCAGCGCCAAGCCAGTTA"]
        spoligo_dictionary["spacer20"] = ["TTGACCTCGCCAGGAGAGAAGATCA", "TGATCTTCTCTCCTGGCGAGGTCAA"]
        spoligo_dictionary["spacer21"] = ["TCGATGTCGATGTCCCAATCGTCGA", "TCGACGATTGGGACATCGACATCGA"]
        spoligo_dictionary["spacer22"] = ["ACCGCAGACGGCACGATTGAGACAA", "TTGTCTCAATCGTGCCGTCTGCGGT"]
        spoligo_dictionary["spacer23"] = ["AGCATCGCTGATGCGGTCCAGCTCG", "CGAGCTGGACCGCATCAGCGATGCT"]
        spoligo_dictionary["spacer24"] = ["CCGCCTGCTGGGTGAGACGTGCTCG", "CGAGCACGTCTCACCCAGCAGGCGG"]
        spoligo_dictionary["spacer25"] = ["GATCAGCGACCACCGCACCCTGTCA", "TGACAGGGTGCGGTGGTCGCTGATC"]
        spoligo_dictionary["spacer26"] = ["CTTCAGCACCACCATCATCCGGCGC", "GCGCCGGATGATGGTGGTGCTGAAG"]
        spoligo_dictionary["spacer27"] = ["GGATTCGTGATCTCTTCCCGCGGAT", "ATCCGCGGGAAGAGATCACGAATCC"]
        spoligo_dictionary["spacer28"] = ["TGCCCCGGCGTTTAGCGATCACAAC", "GTTGTGATCGCTAAACGCCGGGGCA"]
        spoligo_dictionary["spacer29"] = ["AAATACAGGCTCCACGACACGACCA", "TGGTCGTGTCGTGGAGCCTGTATTT"]
        spoligo_dictionary["spacer30"] = ["GGTTGCCCCGCGCCCTTTTCCAGCC", "GGCTGGAAAAGGGCGCGGGGCAACC"]
        spoligo_dictionary["spacer31"] = ["TCAGACAGGTTCGCGTCGATCAAGT", "ACTTGATCGACGCGAACCTGTCTGA"]
        spoligo_dictionary["spacer32"] = ["GACCAAATAGGTATCGGCGTGTTCA", "TGAACACGCCGATACCTATTTGGTC"]
        spoligo_dictionary["spacer33"] = ["GACATGACGGCGGTGCCGCACTTGA", "TCAAGTGCGGCACCGCCGTCATGTC"]
        spoligo_dictionary["spacer34"] = ["AAGTCACCTCGCCCACACCGTCGAA", "TTCGACGGTGTGGGCGAGGTGACTT"]
        spoligo_dictionary["spacer35"] = ["TCCGTACGCTCGAAACGCTTCCAAC", "GTTGGAAGCGTTTCGAGCGTACGGA"]
        spoligo_dictionary["spacer36"] = ["CGAAATCCAGCACCACATCCGCAGC", "GCTGCGGATGTGGTGCTGGATTTCG"]
        spoligo_dictionary["spacer37"] = ["CGCGAACTCGTCCACAGTCCCCCTT", "AAGGGGGACTGTGGACGAGTTCGCG"]
        spoligo_dictionary["spacer38"] = ["CGTGGATGGCGGATGCGTTGTGCGC", "GCGCACAACGCATCCGCCATCCACG"]
        spoligo_dictionary["spacer39"] = ["GACGATGGCCAGTAAATCGGCGTGG", "CCACGCCGATTTACTGGCCATCGTC"]
        spoligo_dictionary["spacer40"] = ["CGCCATCTGTGCCTCATACAGGTCC", "GGACCTGTATGAGGCACAGATGGCG"]
        spoligo_dictionary["spacer41"] = ["GGAGCTTTCCGGCTTCTATCAGGTA", "TACCTGATAGAAGCCGGAAAGCTCC"]
        spoligo_dictionary["spacer42"] = ["ATGGTGGGACATGGACGAGCGCGAC", "GTCGCGCTCGTCCATGTCCCACCAT"]
        spoligo_dictionary["spacer43"] = ["CGCAGAATCGCACCGGGTGCGGGAG", "CTCCCGCACCCGGTGCGATTCTGCG"]

#        spoligo_dictionary["spacer01"] = ["TGATCCAGAGCCGGCGACCCTCTAT"]
#        spoligo_dictionary["spacer02"] = ["CAAAAGCTGTCGCCCAAGCATGAGG"]
#        spoligo_dictionary["spacer03"] = ["CCGTGCTTCCAGTGATCGCCTTCTA"]
#        spoligo_dictionary["spacer04"] = ["ACGTCATACGCCGACCAATCATCAG"]
#        spoligo_dictionary["spacer05"] = ["TTTTCTGACCACTTGTGCGGGATTA"]
#        spoligo_dictionary["spacer06"] = ["CGTCGTCATTTCCGGCTTCAATTTC"]
#        spoligo_dictionary["spacer07"] = ["GAGGAGAGCGAGTACTCGGGGCTGC"]
#        spoligo_dictionary["spacer08"] = ["CGTGAAACCGCCCCCAGCCTCGCCG"]
#        spoligo_dictionary["spacer09"] = ["ACTCGGAATCCCATGTGCTGACAGC"]
#        spoligo_dictionary["spacer10"] = ["TCGACACCCGCTCTAGTTGACTTCC"]
#        spoligo_dictionary["spacer11"] = ["GTGAGCAACGGCGGCGGCAACCTGG"]
#        spoligo_dictionary["spacer12"] = ["ATATCTGCTGCCCGCCCGGGGAGAT"]
#        spoligo_dictionary["spacer13"] = ["GACCATCATTGCCATTCCCTCTCCC"]
#        spoligo_dictionary["spacer14"] = ["GGTGTGATGCGGATGGTCGGCTCGG"]
#        spoligo_dictionary["spacer15"] = ["CTTGAATAACGCGCAGTGAATTTCG"]
#        spoligo_dictionary["spacer16"] = ["CGAGTTCCCGTCAGCGTCGTAAATC"]
#        spoligo_dictionary["spacer17"] = ["GCGCCGGCCCGCGCGGATGACTCCG"]
#        spoligo_dictionary["spacer18"] = ["CATGGACCCGGGCGAGCTGCAGATG"]
#        spoligo_dictionary["spacer19"] = ["TAACTGGCTTGGCGCTGATCCTGGT"]
#        spoligo_dictionary["spacer20"] = ["TTGACCTCGCCAGGAGAGAAGATCA"]
#        spoligo_dictionary["spacer21"] = ["TCGATGTCGATGTCCCAATCGTCGA"]
#        spoligo_dictionary["spacer22"] = ["ACCGCAGACGGCACGATTGAGACAA"]
#        spoligo_dictionary["spacer23"] = ["AGCATCGCTGATGCGGTCCAGCTCG"]
#        spoligo_dictionary["spacer24"] = ["CCGCCTGCTGGGTGAGACGTGCTCG"]
#        spoligo_dictionary["spacer25"] = ["GATCAGCGACCACCGCACCCTGTCA"]
#        spoligo_dictionary["spacer26"] = ["CTTCAGCACCACCATCATCCGGCGC"]
#        spoligo_dictionary["spacer27"] = ["GGATTCGTGATCTCTTCCCGCGGAT"]
#        spoligo_dictionary["spacer28"] = ["TGCCCCGGCGTTTAGCGATCACAAC"]
#        spoligo_dictionary["spacer29"] = ["AAATACAGGCTCCACGACACGACCA"]
#        spoligo_dictionary["spacer30"] = ["GGTTGCCCCGCGCCCTTTTCCAGCC"]
#        spoligo_dictionary["spacer31"] = ["TCAGACAGGTTCGCGTCGATCAAGT"]
#        spoligo_dictionary["spacer32"] = ["GACCAAATAGGTATCGGCGTGTTCA"]
#        spoligo_dictionary["spacer33"] = ["GACATGACGGCGGTGCCGCACTTGA"]
#        spoligo_dictionary["spacer34"] = ["AAGTCACCTCGCCCACACCGTCGAA"]
#        spoligo_dictionary["spacer35"] = ["TCCGTACGCTCGAAACGCTTCCAAC"]
#        spoligo_dictionary["spacer36"] = ["CGAAATCCAGCACCACATCCGCAGC"]
#        spoligo_dictionary["spacer37"] = ["CGCGAACTCGTCCACAGTCCCCCTT"]
#        spoligo_dictionary["spacer38"] = ["CGTGGATGGCGGATGCGTTGTGCGC"]
#        spoligo_dictionary["spacer39"] = ["GACGATGGCCAGTAAATCGGCGTGG"]
#        spoligo_dictionary["spacer40"] = ["CGCCATCTGTGCCTCATACAGGTCC"]
#        spoligo_dictionary["spacer41"] = ["GGAGCTTTCCGGCTTCTATCAGGTA"]
#        spoligo_dictionary["spacer42"] = ["ATGGTGGGACATGGACGAGCGCGAC"]
#        spoligo_dictionary["spacer43"] = ["CGCAGAATCGCACCGGGTGCGGGAG"]

        count_summary={}

        with futures.ProcessPoolExecutor() as pool: #max_workers=4
            for v, count in pool.map(script1.finding_sp, spoligo_dictionary.values()):
                for k, value in spoligo_dictionary.items():
                    if v == value:
                        count_summary.update({k:count})
                        count_summary=OrderedDict(sorted(count_summary.items()))

#        for v in spoligo_dictionary.values():
#            v, count = script1.finding_sp(v)
#            print(v)
#            print(total)
#            for k, value in spoligo_dictionary.items():
#                if v == value:
#                    count_summary.update({k:count})
#                    count_summary=OrderedDict(sorted(count_summary.items()))

        #Binary dictionary
        spoligo_binary_dictionary={}
        for k, v in count_summary.items():
            if v > 4:
                spoligo_binary_dictionary.update({k:1})
            else:
                spoligo_binary_dictionary.update({k:0})
        spoligo_binary_dictionary=OrderedDict(sorted(spoligo_binary_dictionary.items()))
        
        spoligo_binary_list=[]
        for v in spoligo_binary_dictionary.values():
            spoligo_binary_list.append(v)
        bovis_string=''.join(str(e) for e in spoligo_binary_list) #bovis_string correct
        hexadecimal = script1.binary_to_hex(bovis_string)
        
        write_out = open("spoligo.txt", 'w')
        
        found = False
        with open(spoligo_db) as f: # put into dictionary or list
            for line in f:
                line=line.rstrip()
                octalcode = line.split()[0] #no arg splits on whitespace
                sbcode = line.split()[1]
                binarycode = line.split()[2]
                if bovis_string == '0000000000000000000000000000000000000000000':
                    found=True
                    print("%s %s %s %s" % (octalcode, sbcode, hexadecimal, binarycode))
                    print("%s %s %s %s" % (octalcode, sbcode, hexadecimal, binarycode), file=write_out)
                    print("CHECK SAMPLE!  NO SPACERS FOUND.  LIKELY NOT TB COMPLEX.  ALTHOUGH SB2277 IS A ZERO STRING BINARY\n")
                    print("CHECK SAMPLE!  NO SPACERS FOUND.  LIKELY NOT TB COMPLEX.  ALTHOUGH SB2277 IS A ZERO STRING BINARY", file=write_out)
                    print("\nOne mismatch allowed spacers search against both R1 and R2 reads.\n", file=write_out)
                    for k, v in count_summary.items():
                        print(k, v, file=write_out)
                elif bovis_string == binarycode:
                    found=True
                    print("Pattern found:")
                    print("%s %s %s %s" % (octalcode, sbcode, hexadecimal, binarycode))
                    print("%s %s %s %s" % (octalcode, sbcode, hexadecimal, binarycode), file=write_out)
                    print("\One mismatch allowed spacer search against both R1 and R2 reads.\n", file=write_out)
                    for k, v in count_summary.items():
                        print(k, v, file=write_out)
        if not found:
            octal = script1.binary_to_octal(bovis_string)
            sbcode = "N/A"
            print("%s %s %s %s" % (octal, sbcode, hexadecimal, bovis_string))
            print("%s %s %s %s" % (octal, sbcode, hexadecimal, bovis_string), file=write_out)
            print("SPOLIGO SB NUMBER NOT FOUND\n")
            print("\nSPOLIGO SB NUMBER NOT FOUND\n", file=write_out)
            print("\nOne mismatch allowed spacer search against both R1 and R2 reads.\n", file=write_out)
            
            for k, v in count_summary.items():
                print(k, v, file=write_out)

        print("bovis_string: %s" % bovis_string, file=write_out)
        print("binarycode  : %s" % binarycode, file=write_out)

        for i in fastqs: #remove unzipped fastq files to save space
            os.remove(i)
            
        write_out.close()
###############################################

    def best_reference(self):
    
        '''Use oligos to determine species.  Most often if the absents of a single oligo from a set specific for either brucella or bovis will confer species type.  Some species will the absents of more than one oligo.  Oligo findings are translated to binary patterns.'''
        
        fastqs = glob.glob(zips + '/*.fastq')
        if len(fastqs) < 2:
            script1.unzipfiles()
        fastqs = glob.glob(zips + '/*.fastq')

        print("\nFinding the best reference\n")
        
        write_out = open("best_reference.txt", 'w')
      
        '''get the species'''
        oligo_dictionary = {}
        oligo_dictionary["01_ab1"] = "AATTGTCGGATAGCCTGGCGATAACGACGC"
        oligo_dictionary["02_ab3"] = "CACACGCGGGCCGGAACTGCCGCAAATGAC"
        oligo_dictionary["03_ab5"] = "GCTGAAGCGGCAGACCGGCAGAACGAATAT"
        oligo_dictionary["04_mel"] = "TGTCGCGCGTCAAGCGGCGTGAAATCTCTG"
        oligo_dictionary["05_suis1"] = "TGCGTTGCCGTGAAGCTTAATTCGGCTGAT"
        oligo_dictionary["06_suis2"] = "GGCAATCATGCGCAGGGCTTTGCATTCGTC"
        oligo_dictionary["07_suis3"] = "CAAGGCAGATGCACATAATCCGGCGACCCG"
        oligo_dictionary["08_ceti1"] = "GTGAATATAGGGTGAATTGATCTTCAGCCG"
        oligo_dictionary["09_ceti2"] = "TTACAAGCAGGCCTATGAGCGCGGCGTGAA"
        oligo_dictionary["10_canis4"] = "CTGCTACATAAAGCACCCGGCGACCGAGTT"
        oligo_dictionary["11_canis"] = "ATCGTTTTGCGGCATATCGCTGACCACAGC"
        oligo_dictionary["12_ovis"] = "CACTCAATCTTCTCTACGGGCGTGGTATCC"
        oligo_dictionary["13_ether2"] = "CGAAATCGTGGTGAAGGACGGGACCGAACC"
        oligo_dictionary["14_63B1"] = "CCTGTTTAAAAGAATCGTCGGAACCGCTCT"
        oligo_dictionary["15_16M0"] = "TCCCGCCGCCATGCCGCCGAAAGTCGCCGT"
        oligo_dictionary["16_tb157"] = "CTCTTCGTATACCGTTCCGTCGTCACCATGGTCCT"
        oligo_dictionary["17_tb7"] = "TCACGCAGCCAACGATATTCGTGTACCGCGACGGT"
        oligo_dictionary["18_tbbov"] = "CTGGGCGACCCGGCCGACCTGCACACCGCGCATCA"
        oligo_dictionary["19_tb5"] = "CCGTGGTGGCGTATCGGGCCCCTGGATCGCGCCCT"
        oligo_dictionary["20_tb2"] = "ATGTCTGCGTAAAGAAGTTCCATGTCCGGGAAGTA"
        oligo_dictionary["21_tb3"] = "GAAGACCTTGATGCCGATCTGGGTGTCGATCTTGA"
        oligo_dictionary["22_tb4"] = "CGGTGTTGAAGGGTCCCCCGTTCCAGAAGCCGGTG"
        oligo_dictionary["23_tb6"] = "ACGGTGATTCGGGTGGTCGACACCGATGGTTCAGA"
        oligo_dictionary["24_para"] = "CCTTTCTTGAAGGGTGTTCG"
        oligo_dictionary["25_para2"] = "CGAACACCCTTCAAGAAAGG"

        brucella_identifications = {}
        brucella_identifications["111111111111111"] = "odd" #Unexpected findings
        brucella_identifications["011111111111111"] = "ab1" #Brucella abortus bv 1, 2 or 4
        brucella_identifications["101111111111111"] = "ab3" #Brucella abortus bv 3
        brucella_identifications["110111111111111"] = "ab1" #Brucella abortus bv 5, 6 or 9
        brucella_identifications["111011111111110"] = "mel1"
        brucella_identifications["000001010110110"] = "mel1"
        brucella_identifications["111011111111101"] = "mel2"
        brucella_identifications["000001010110100"] = "mel2"
        brucella_identifications["010001010110100"] = "mel2"
        brucella_identifications["111001111110101"] = "mel2"
        brucella_identifications["111011111111011"] = "mel3"
        brucella_identifications["111001111110011"] = "mel3"
        brucella_identifications["111101111111111"] = "suis1"
        brucella_identifications["111110111111111"] = "suis2"
        brucella_identifications["111111011111110"] = "suis3"
        brucella_identifications["111111101111111"] = "ceti1"
        brucella_identifications["111111100111111"] = "ceti1"
        brucella_identifications["111111110111111"] = "ceti2"
        brucella_identifications["111111111011110"] = "suis4"
        brucella_identifications["111111111001110"] = "canis"
        brucella_identifications["111111111110111"] = "ovis"

        bovis_identifications = {}
        bovis_identifications["11101111"] = "h37" #tb1
        bovis_identifications["11101101"] = "h37" #tb1
        bovis_identifications["01100111"] = "h37" #tb2
        bovis_identifications["01101011"] = "h37" #tb3
        bovis_identifications["11101011"] = "h37" #tb3
        bovis_identifications["01101111"] = "h37" #tb4a
        bovis_identifications["01101101"] = "h37" #tb4b
        bovis_identifications["11101101"] = "h37" #tb4b
        bovis_identifications["01101111"] = "h37" #tb4b
        bovis_identifications["11111111"] = "h37" #tb5
        bovis_identifications["11001111"] = "h37" #tb6
        bovis_identifications["10101110"] = "h37" #tb7
        bovis_identifications["11001110"] = "bovis" #bovis
        bovis_identifications["11011110"] = "bovis" #bovis
        bovis_identifications["11001100"] = "bovis" #bovis
        
        para_identifications = {}
        para_identifications["1"] = "para"
        para_identifications["01"] = "para"
        para_identifications["11"] = "para"

        count_summary={}

        with futures.ProcessPoolExecutor() as pool: #max_workers=4
            for v, count in pool.map(script1.finding_best_ref, oligo_dictionary.values()):
                for k, value in oligo_dictionary.items():
                    if v == value:
                        count_summary.update({k:count})
                        count_summary=OrderedDict(sorted(count_summary.items()))

        count_list=[]
        for v in count_summary.values():
            count_list.append(v)
        brucella_sum=sum(count_list[:15])
        bovis_sum=sum(count_list[15:23])
        para_sum=sum(count_list[23:])
        
        print("Best reference Brucella counts:", file=write_out)
        for i in count_list[:15]:
            print(i,  end=',', file=write_out)
            
        print("\nBest reference TB counts:", file=write_out)
        for i in count_list[15:23]:
            print(i,  end=',', file=write_out)

        print("\nBest reference Para counts:", file=write_out)
        for i in count_list[23:]:
            print(i,  end=',', file=write_out)

        #Binary dictionary
        binary_dictionary={}
        for k, v in count_summary.items():
            if v > 1:
                binary_dictionary.update({k:1})
            else:
                binary_dictionary.update({k:0})
        binary_dictionary=OrderedDict(sorted(binary_dictionary.items()))

        binary_list=[]
        for v in binary_dictionary.values():
            binary_list.append(v)
        brucella_binary=binary_list[:15]
        brucella_string=''.join(str(e) for e in brucella_binary)
        bovis_binary=binary_list[15:23]
        bovis_string=''.join(str(e) for e in bovis_binary)
        para_binary=binary_list[23:]
        para_string=''.join(str(e) for e in para_binary)

        if brucella_sum > 3:
            if brucella_string in brucella_identifications:
                print("Brucella group, species %s" % brucella_identifications[brucella_string])
                print("\n\nBrucella group, species %s" % brucella_identifications[brucella_string], file=write_out)
                return(brucella_identifications[brucella_string]) # return to set parameters
            else:
                print("Brucella group, but no match")
                print("\n\nBrucella group, but no match", file=write_out)
        elif bovis_sum > 3:
            if bovis_string in bovis_identifications:
                print("TB group, species %s" % bovis_identifications[bovis_string])
                print("\n\nTB group, species %s" % bovis_identifications[bovis_string], file=write_out)
                return(bovis_identifications[bovis_string]) # return to set parameters
            else:
                print("TB group, but no match")
                print("\n\nTB group, but no match", file=write_out)
        elif para_sum >= 1:
            if para_string in para_identifications:
                print("Para group")
                print("\n\nPara group", file=write_out)
                return("para") # return to set parameters
            else:
                print("No match")
                print("\n\nNo match", file=write_out)

        write_out.close()
        
        for i in fastqs: #remove unzipped fastq files to save space
            os.remove(i)
###############################################

    def add_zero_coverage(coverage_in, vcf_file, loc_sam):
        
        temp_vcf = loc_sam + "-temp.vcf"
        zero_coverage_vcf = loc_sam + "_zc.vcf"
        
        zero_position=[]
        total_length = 0
        total_zero_coverage = 0
        with open(coverage_in) as f:
            for line in f:
                total_length = total_length + 1
                line.rstrip()
                line=re.split(':|\t', line)
                chromosome=line[0]
                position=line[1]
                abs_pos = chromosome + "-" + position
                depth=line[2]
                if depth == "0":
                    zero_position.append(abs_pos) #positions with zero coverage in coverage file
                    total_zero_coverage = total_zero_coverage + 1
            print(len(zero_position))
        
        genome_coverage = 0
        total_coverage = total_length - total_zero_coverage

        genome_coverage =  "{:.2%}".format(total_coverage/total_length)

        average_list = []
        with open(coverage_in) as f:
            for line in f:
                line.rstrip()
                line=re.split(':|\t', line)
                depth=str(line[2])
                if depth.isdigit():
                    depth = int(depth)
                    average_list.append(depth)
            ave_coverage = mean(average_list)

        zero_position_found=[]
        write_out=open(temp_vcf, 'w')
        with open(vcf_file) as f:
            for line in f:
                line=line.rstrip()
                if line.startswith("#"): # save headers to file
                    print(line, file=write_out)
                elif not line.startswith("#"): # position rows
                    split_line = line.split('\t')
                    chromosome=split_line[0] # get chromosome
                    position=split_line[1] # get position
                    abs_pos = chromosome + "-" + position
                    ref=split_line[3] # REF
                    alt=split_line[4] # ALT
                    ref_len=len(ref)
                    alt_len=len(alt)
                    if abs_pos in zero_position: # if a position has zero coverage
                        print("%s is in zeropostions" % position)
                        zero_position_found.append(position)
                        print("%s\t%s\t.\tN\t.\t.\t.\t.\tGT\t./." % (chromosome, position), file=write_out) # print a zero coverage line
                    elif ref_len == 1 and alt_len == 1:
                        print(line, file=write_out)
            print("##### Chromosome: %s" % chromosome)
            #zero_position = list(map(int, zero_position)) # change list items from string to numbers
            #zero_position_found = list(map(int, zero_position_found))

            print("using list comprehension")
            zero_not_added = [x for x in zero_position if x not in  zero_position_found] # use list comprehension to subtract one list from the other
            for abs_position in zero_not_added:
                split_line = abs_position.split('-')
                chromosome=split_line[0]
                position=split_line[1]
                print("%s\t%s\t.\tN\t.\t.\t.\t.\tGT\t./." % (chromosome, position), file=write_out) # print a zero coverage line
        write_out.close()

        os.system("picard SortVcf INPUT={} OUTPUT={}" .format(temp_vcf, zero_coverage_vcf))
        #os.system("vcf-sort {} > {}" .format(temp_vcf, zero_coverage_vcf))
        os.remove(temp_vcf)
        
        vcf_reader = vcf.Reader(open(zero_coverage_vcf, 'r'))
        good_snp_count = 0
        for record in vcf_reader:
            try:
                chrom = record.CHROM
                position = record.POS
                try:
                    if str(record.ALT[0]) != "None" and record.INFO['AC'][0] == 2 and len(record.REF) == 1 and record.QUAL > 150:
                        good_snp_count = good_snp_count + 1
                except KeyError:
                    pass
            except ZeroDivisionError:
                print ("ZeroDivisionError error found")
            except ValueError:
                print ("ValueError error found")
            except UnboundLocalError:
                print ("UnboundLocalError error found")
            except TypeError:
                print ("TypeError error found")

        return(zero_coverage_vcf, good_snp_count, ave_coverage, genome_coverage)

###############################################

    def align_reads(self):
    
        if self.species == "NO FINDINGS":
            read_base = os.path.basename(R1)
            sample_name=re.sub('_.*', '', read_base)
            R1size = script1.sizeof_fmt(os.path.getsize(R1))
            R2size = script1.sizeof_fmt(os.path.getsize(R2))
            stat_summary = {}
            stat_summary["01"] = sample_name
            stat_summary["02"] = "NOT_FOUND"
            stat_summary["03"] = "N/A"
            stat_summary["04"] = R1size
            stat_summary["05"] = R2size
            stat_summary["06"] = "CHECK SAMPLE *****************************************"
            stat_summary=OrderedDict(sorted(stat_summary.items()))
            return(stat_summary)
        else:
            startTime = datetime.now()
            print ("\n\n*** START ***\n")
            print ("Start time: %s" % startTime)

            read_base = os.path.basename(R1)
            sample_name=re.sub('_.*', '', read_base)

            sample_reference=glob.glob(directory + '/*fasta')
            hqs=glob.glob(directory + '/*vcf')
            
            print("self.species: %s" % self.species)
            if self.species in ["ab1", "ab3", "suis1", "suis3", "suis4", "mel1", "mel2", "mel3", "canis", "ceti1", "ceti2"]:
                print("Brucella")
                self.mlst()
            elif self.species in ["bovis", "h37", "af"]:
                print("TB")
                self.spoligo()
            
            print ("reference: %s" % sample_reference)
            ref=re.sub('\.fasta', '', os.path.basename(sample_reference[0]))
            if len(sample_reference) != 1:
                print("### ERROR reference not available or too many")
                sys.exit(0)
            if len(hqs) != 1:
                print("### ERROR high quality snps not available or too many")
                sys.exit(0)
            sample_reference=sample_reference[0]
            hqs=hqs[0]
            
            print ("--")
            print("sample name: %s" % sample_name)
            print("sample reference: %s" % sample_reference)
            print("Read 1: %s" % R1)
            print("Read 2: %s\n" % R2)
            print("directory: %s" % directory)
            print ("--")

            loc_sam=directory + "/" + sample_name
            
            os.system("samtools faidx {}" .format(sample_reference))
            os.system("picard CreateSequenceDictionary REFERENCE={} OUTPUT={}" .format(sample_reference, directory + "/" + ref + ".dict"))
            os.system("bwa index {}" .format(sample_reference))
            samfile = loc_sam + ".sam"
            allbam = loc_sam + "-all.bam"
            unmapsam = loc_sam + "-unmapped.sam"
            unmapped_read1 = loc_sam + "-unmapped_R1.fastq"
            unmapped_read2 = loc_sam + "-unmapped_R2.fastq"
            unmapped_read1gz = loc_sam + "-unmapped_R1.fastq.gz"
            unmapped_read2gz = loc_sam + "-unmapped_R2.fastq.gz"
            abyss_out = loc_sam + "-unmapped_contigs.fasta"
            sortedbam = loc_sam + "-sorted.bam"
            nodupbam = loc_sam + "-nodup.bam"
            metrics = loc_sam + "-metrics.txt"
            indel_realigner = loc_sam + ".intervals"
            realignedbam = loc_sam + "-realigned.bam"
            recal_group = loc_sam + "-recal_group"
            prebam=loc_sam + "-pre.bam"
            qualitybam = loc_sam + "-quality.bam"
            coverage_file=loc_sam + "-coverage.txt"
            hapall = loc_sam + "-hapall.vcf"
            bamout = loc_sam + "-bamout.bam"
            
            print("\n@@@ BWA mem")
            os.system("bwa mem -M -t 16 {} {} {} > {}" .format(sample_reference, R1, R2, samfile))

            print("\nAdd read group and out all BAM")
            os.system("picard AddOrReplaceReadGroups INPUT={} OUTPUT={} RGLB=lib1 RGPU=unit1 RGSM={} RGPL=illumina" .format(samfile, allbam, sample_name))
            os.system("samtools index {}" .format(allbam))

            print("\n@@@ Samtools unmapped")
            os.system("samtools view -h -f4 -T {} {} -o {}" .format(sample_reference, allbam, unmapsam))

            print("\n@@@ Unmapped to FASTQ")
            os.system("picard SamToFastq INPUT={} FASTQ={} SECOND_END_FASTQ={}" .format(unmapsam, unmapped_read1, unmapped_read2))
            
            print("\n@@@ Abyss")
            abyss_contig_count=0
<<<<<<< HEAD
            os.system("ABYSS --out {} --coverage 3 --kmer 64 {} {}" .format(abyss_out, unmapped_read1, unmapped_read2))
            with open(abyss_out) as f:
                for line in f:
                    abyss_contig_count += line.count(">")
=======
            os.system("ABYSS --out {} --kmer 64 {} {}" .format(abyss_out, unmapped_read1, unmapped_read2))
            try:
                with open(abyss_out) as f:
                    for line in f:
                        abyss_contig_count += line.count(">")
            except FileNotFoundError:
                abyss_contig_count = 0
>>>>>>> origin/master

            print("\n@@@ Sort BAM")
            os.system("samtools sort {} -o {}" .format(allbam, sortedbam))
            os.system("samtools index {}" .format(sortedbam))
            
            print("\n@@@ Write stats to file")
            stat_file = "stat_align.txt"
            stat_out = open(stat_file, 'w')
            #os.system("samtools idxstats {} > {}" .format(sortedbam, stat_out)) Doesn't work when needing to std out.
            stat_out.write(os.popen("samtools idxstats {} " .format(sortedbam)).read())
            stat_out.close()

            with open(stat_file) as f:
                first_line = f.readline()
                first_line = first_line.rstrip()
                first_line=re.split(':|\t', first_line)
                reference_sequence_name = str(first_line[0])
                sequence_length = "{:,}".format(int(first_line[1]))
                allbam_mapped_reads = int(first_line[2])
                allbam_unmapped_reads = "{:,}".format(int(first_line[3]))

            print("\n@@@ Find duplicate reads")
            os.system("picard MarkDuplicates INPUT={} OUTPUT={} METRICS_FILE={} ASSUME_SORTED=true REMOVE_DUPLICATES=true" .format(sortedbam, nodupbam, metrics))
            os.system("samtools index {}" .format(nodupbam))
            
            duplicate_stat_file = "duplicate_stat_align.txt"
            duplicate_stat_out = open(duplicate_stat_file, 'w')
            #os.system("samtools idxstats {} > {}" .format(sortedbam, stat_out)) Doesn't work when needing to std out.
            duplicate_stat_out.write(os.popen("samtools idxstats {} " .format(nodupbam)).read())
            duplicate_stat_out.close()
            with open(duplicate_stat_file) as f:
                dup_first_line = f.readline()
                dup_first_line = dup_first_line.rstrip()
                dup_first_line=re.split(':|\t', dup_first_line)
                nodupbam_mapped_reads = int(dup_first_line[2])
                nodupbam_unmapped_reads = int(dup_first_line[3])
            try:
                unmapped_reads = allbam_mapped_reads - nodupbam_mapped_reads
            except:
                unmapped_reads = "none_found"
            
            allbam_mapped_reads = "{:,}".format(allbam_mapped_reads)
            print(unmapped_reads)

            print("\n@@@  Realign indels")
            os.system("gatk -T RealignerTargetCreator -I {} -R {} -o {}" .format(nodupbam, sample_reference, indel_realigner))
            if not os.path.isfile(indel_realigner):
                os.system("gatk -T RealignerTargetCreator --fix_misencoded_quality_scores -I {} -R {} -o {}" .format(nodupbam, sample_reference, indel_realigner))
            os.system("gatk -T IndelRealigner -I {} -R {} -targetIntervals {} -o {}" .format(nodupbam, sample_reference, indel_realigner, realignedbam))
            if not os.path.isfile(realignedbam):
                os.system("gatk -T IndelRealigner --fix_misencoded_quality_scores -I {} -R {} -targetIntervals {} -o {}" .format(nodupbam, sample_reference, indel_realigner, realignedbam))

            print("\n@@@ Base recalibration")
            os.system("gatk -T BaseRecalibrator -I {} -R {} -knownSites {} -o {}". format(realignedbam, sample_reference, hqs, recal_group))
            if not os.path.isfile(realignedbam):
                os.system("gatk -T BaseRecalibrator  --fix_misencoded_quality_scores -I {} -R {} -knownSites {} -o {}". format(realignedbam, sample_reference, hqs, recal_group))

            print("\n@@@ Make realigned BAM")
            os.system("gatk -T PrintReads -R {} -I {} -BQSR {} -o {}" .format (sample_reference, realignedbam, recal_group, prebam))
            if not os.path.isfile(prebam):
                os.system("gatk -T PrintReads  --fix_misencoded_quality_scores -R {} -I {} -BQSR {} -o {}" .format (sample_reference, realignedbam, recal_group, prebam))

            print("\n@@@ Clip reads")
            os.system("gatk -T ClipReads -R {} -I {} -o {} -filterNoBases -dcov 10" .format(sample_reference, prebam, qualitybam))
            os.system("samtools index {}" .format(qualitybam))

            print("\n@@@ Depth of coverage using GATK")
            os.system("gatk -T DepthOfCoverage -R {} -I {} -o {} -omitIntervals --omitLocusTable --omitPerSampleStats -nt 8" .format(sample_reference, prebam, coverage_file))

            print("\n@@@ Calling SNPs with HaplotypeCaller")
            os.system("gatk -R {} -T HaplotypeCaller -I {} -o {} -bamout {} -dontUseSoftClippedBases -allowNonUniqueKmersInRef" .format(sample_reference, qualitybam, hapall, bamout))

            print("Getting Zero Coverage...\n")
            zero_coverage_vcf, good_snp_count, ave_coverage, genome_coverage = script1.add_zero_coverage(coverage_file, hapall, loc_sam)

            ###
            if gbk_file is not "None":
                try:
                    in_annotation_as_dict = SeqIO.to_dict(SeqIO.parse(gbk_file, "genbank"))
                    annotated_vcf = loc_sam + "-annotated.vcf"
                    write_out=open(annotated_vcf, 'w')
                    
                    with open(zero_coverage_vcf) as vfile:
                        print("finding annotations...\n")
                        for line in vfile:
                            annotated_line = script1.get_annotations(line, in_annotation_as_dict)
                            print("%s" % annotated_line, file=write_out)
                    write_out.close()
                except AttributeError:
                    pass
#
#            with open(zero_coverage_vcf) as vfile:
#                print("finding annotations...\n")
#                with futures.ProcessPoolExecutor() as pool:
#                    for annotated_line in pool.map(partial(script1.get_annotations, in_annotation_as_dict = in_annotation_as_dict), vfile): #function, passed arguement, list
#                        print("%s" % annotated_line, file=write_out)
#            write_out.close()
            ###

            os.remove(samfile)
            os.remove(allbam)
            os.remove(nodupbam)
            os.remove(nodupbam + ".bai")
            os.remove(unmapsam)
            os.remove(sortedbam)
            os.remove(sortedbam + ".bai")
            os.remove(indel_realigner)
            os.remove(realignedbam)
            os.remove(loc_sam + "-realigned.bai")
            os.remove(recal_group)
            os.remove(prebam)
            os.remove(loc_sam + "-pre.bai")
            os.remove(hqs)
            os.remove(hqs + ".idx")
            os.remove(sample_reference + ".amb")
            os.remove(sample_reference + ".ann")
            os.remove(sample_reference + ".bwt")
            os.remove(sample_reference + ".pac")
            os.remove(sample_reference + ".sa")
            os.remove(ref + ".dict")
            os.remove(duplicate_stat_file)
            os.remove(stat_file)

            unmapped = directory + "/unmapped"
            os.makedirs(unmapped)

            newZip = zipfile.ZipFile(unmapped_read1gz, 'w')
            newZip.write(unmapped_read1, compress_type=zipfile.ZIP_DEFLATED)
            newZip.close()
            newZip = zipfile.ZipFile(unmapped_read2gz, 'w')
            newZip.write(unmapped_read2, compress_type=zipfile.ZIP_DEFLATED)
            newZip.close()
            os.remove(unmapped_read1)
            os.remove(unmapped_read2)
            
            try:
                shutil.move(unmapped_read1gz, unmapped)
                shutil.move(unmapped_read2gz, unmapped)
                shutil.move(abyss_out, unmapped)
            except FileNotFoundError:
                pass
            
            alignment = directory + "/alignment"
            os.makedirs(alignment)
            movefiles = glob.glob('*-*')
            for i in movefiles:
                shutil.move(i, alignment)
            try:
                shutil.move(sample_reference, alignment)
                shutil.move(sample_reference + ".fai", alignment)
            except shutil.Error:
                pass
            except FileNotFoundError:
                pass
            except FileExistsError:
                pass

            runtime = (datetime.now() - startTime)
            print ("\n\nruntime: %s:  \n" % runtime)
            ave_coverage = "{:0.1f}".format(float(ave_coverage))
            print("average_coverage: %s" % ave_coverage)

            R1size = script1.sizeof_fmt(os.path.getsize(R1))
            R2size = script1.sizeof_fmt(os.path.getsize(R2))

            try:
                with open("mlst.txt") as f:
                    first_line = f.readline()
                    mlst_type = first_line.rstrip()
                    first_line = first_line.split()
                    mlst_type = first_line[1:]
                    mlst_type = '-'.join(mlst_type)
            except FileNotFoundError:
                mlst_type = "N/A"

            try:
                with open("spoligo.txt") as f:
                    first_line = f.readline()
                    first_line = first_line.rstrip()
                    first_line = first_line.split()
                    octalcode = first_line[0]
                    sbcode = first_line[1]
                    hexcode = first_line[2]
                    binarycode = first_line[3]
            except FileNotFoundError:
                octalcode = "N/A"
                sbcode = "N/A"
                hexcode = "N/A"
                binarycode = "N/A"

            sequence_count = 0
            total_length = 0
            with gzip.open(R2, "rt") as handle:
                for r in SeqIO.parse(handle, "fastq"):
                    total_length = total_length + len(r.seq)
                    sequence_count = sequence_count + 1
            ave_read_length = total_length/sequence_count
            ave_read_length = "{:0.1f}".format(float(ave_read_length))

            stat_summary={}
            stat_summary["01-sample_name"] = sample_name
            stat_summary["02-self.species"] = self.species
            stat_summary["03-reference_sequence_name"] = reference_sequence_name
            stat_summary["04-R1size"] = R1size
            stat_summary["05-R2size"] = R2size
            stat_summary["06-allbam_mapped_reads"] = allbam_mapped_reads
            stat_summary["07-genome_coverage"] = genome_coverage
            stat_summary["08-ave_coverage"] = ave_coverage
            stat_summary["09-ave_read_length"] = ave_read_length
            stat_summary["10-unmapped_reads"] = unmapped_reads
            stat_summary["11-abyss_contig_count"] = abyss_contig_count
            stat_summary["12-good_snp_count"] = good_snp_count
            stat_summary["13-mlst_type"] = mlst_type
            stat_summary["14-octalcode"] = octalcode
            stat_summary["15-sbcode"] = sbcode
            stat_summary["16-hexcode"] = hexcode
            stat_summary["17-binarycode"] = binarycode

            stat_summary=OrderedDict(sorted(stat_summary.items()))

            for k, v in stat_summary.items():
                print("%s: %s" % (k, v))
            
            ###
            # Create a sample stats file in the sample's script1 directory
            ts = time.time()
            st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')
            summary_file = loc_sam + "_" + st + '.xlsx'
            workbook = xlsxwriter.Workbook(summary_file)
            worksheet = workbook.add_worksheet()
            row = 0
            col = 0

            top_row_header = ["sample_name", "self.species", "reference_sequence_name", "R1size", "R2size", "allbam_mapped_reads", "genome_coverage", "ave_coverage", "ave_read_length", "unmapped_reads", "unmapped_assembled_contigs", "good_snp_count", "mlst_type", "octalcode", "sbcode", "hexadecimal_code", "binarycode"]
            for header in top_row_header:
                worksheet.write(row, col, header)
                col += 1
            col = 0
            row += 1
            for v in stat_summary.values():
                    worksheet.write(row, col, v)
                    col += 1
            workbook.close()
            ###
            
            
            return(stat_summary)
###############################################
#R1 = glob.glob('*_R1*fastq.gz')
#R2 = glob.glob('*_R2*fastq.gz')
#sample = script1(R1[0], R2[0], "ab1") #no species give, will find best reference
###sample.show_fastqs()
####sample.best_reference()
####sample.mlst()
####sample.spoligo()
#sample.align_reads()

'''
add para to script 1 and 2

'''

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from concurrent import futures
import multiprocessing
import regex
from collections import OrderedDict
from functools import partial

def finding_sp(R1unzip, R2unzip, v):
    total=0
    total_finds=0
    #if total < 6: # doesn't make a big different.  Might as well get full counts
    # total += sum(seq.count(x) for x in (v)) #v=list of for and rev spacer
    total_finds = [len(regex.findall("(" + spacer + "){s<=1}", seq_string)) for spacer in v]
    for number in total_finds:
        total += number
    return(v, total)

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

def spoligo(R1unzip, R2unzip):

    print("\nFinding spoligotype pattern...\n")
    
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

    count_summary={}

    global seq_string
    sequence_list = []
    for fastq in R1unzip, R2unzip:
        with open(fastq) as in_handle:
            # all 3, title and seq and qual, were needed
            for title, seq, qual in FastqGeneralIterator(in_handle):
                sequence_list.append(seq)
    seq_string = "".join(sequence_list)

#    with futures.ProcessPoolExecutor() as pool:
#        for v, count in pool.map(finding_sp, spoligo_dictionary.values()):
#            for k, value in spoligo_dictionary.items():
#                if v == value:
#                    count_summary.update({k:count})
#                    count_summary=OrderedDict(sorted(count_summary.items()))

    pool = multiprocessing.Pool()
    func = partial(finding_sp, R1unzip, R2unzip)
    for v, count in pool.map(func, spoligo_dictionary.values()):
        for k, value in spoligo_dictionary.items():
            if v == value:
                count_summary.update({k:count})
                count_summary=OrderedDict(sorted(count_summary.items()))
    pool.close()
    pool.join()

#    pool = multiprocessing.Pool()
#    func = partial(finding_best_ref, R1unzip, R2unzip)
#    for v, count in pool.map(func, oligo_dictionary.values()):
#        for k, value in oligo_dictionary.items():
#            if v == value:
#                count_summary.update({k:count})
#                count_summary=OrderedDict(sorted(count_summary.items()))
#    pool.close()
#    pool.join()


    seq_string = ""

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
    hexadecimal = binary_to_hex(bovis_string)
    
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
                octalcode = "spoligo not found"
                sbcode = "spoligo not found"
                hexadecimal = "SB2277 ???"
                binarycode = "0000000000000000000000000000000000000000000"
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
        octal = binary_to_octal(bovis_string)
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

########
########

def finding_best_ref(R1unzip, R2unzip, v):
    count=0
    for fastq in R1unzip, R2unzip:
        with open(fastq) as in_handle:
            # all 3, title and seq and qual, were needed
            for title, seq, qual in FastqGeneralIterator(in_handle):
                count += seq.count(v)
    return(v, count)

def best_reference(R1unzip, R2unzip):
    
#    global R1unzip
#    global R2unzip
#    R1unzip = R1un
#    R2unzip = R2un

    '''Use oligos to determine species.  Most often if the absents of a single oligo from a set specific for either brucella or bovis will confer species type.  Some species will the absents of more than one oligo.  Oligo findings are translated to binary patterns.'''

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
    oligo_dictionary["16_mel1b"] = "TCTGTCCAAACCCCGTGACCGAACAATAGA" #added 2018-01-30
    oligo_dictionary["17_tb157"] = "CTCTTCGTATACCGTTCCGTCGTCACCATGGTCCT"
    oligo_dictionary["18_tb7"] = "TCACGCAGCCAACGATATTCGTGTACCGCGACGGT"
    oligo_dictionary["19_tbbov"] = "CTGGGCGACCCGGCCGACCTGCACACCGCGCATCA"
    oligo_dictionary["20_tb5"] = "CCGTGGTGGCGTATCGGGCCCCTGGATCGCGCCCT"
    oligo_dictionary["21_tb2"] = "ATGTCTGCGTAAAGAAGTTCCATGTCCGGGAAGTA"
    oligo_dictionary["22_tb3"] = "GAAGACCTTGATGCCGATCTGGGTGTCGATCTTGA"
    oligo_dictionary["23_tb4"] = "CGGTGTTGAAGGGTCCCCCGTTCCAGAAGCCGGTG"
    oligo_dictionary["24_tb6"] = "ACGGTGATTCGGGTGGTCGACACCGATGGTTCAGA"
    oligo_dictionary["25_para"] = "CCTTTCTTGAAGGGTGTTCG"
    oligo_dictionary["26_para2"] = "CGAACACCCTTCAAGAAAGG"

    brucella_identifications = {}
    brucella_identifications["1111111111111111"] = "odd" #Unexpected findings
    brucella_identifications["0111111111111111"] = "ab1" #Brucella abortus bv 1, 2 or 4
    brucella_identifications["1011111111111111"] = "ab3" #Brucella abortus bv 3
    brucella_identifications["1101111111111111"] = "ab1" #Brucella abortus bv 5, 6 or 9
    brucella_identifications["1110111111111101"] = "mel1"
    brucella_identifications["0000010101101101"] = "mel1"
    brucella_identifications["1110111111111100"] = "mel1b" #added 2018-01-30
    brucella_identifications["0000010101101100"] = "mel1b" #added 2018-01-30
    brucella_identifications["1110111111111011"] = "mel2"
    brucella_identifications["0000010101101001"] = "mel2"
    brucella_identifications["0100010101101001"] = "mel2"
    brucella_identifications["1110011111101011"] = "mel2"
    brucella_identifications["1110111111110111"] = "mel3"
    brucella_identifications["1110011111100111"] = "mel3"
    brucella_identifications["1111011111111111"] = "suis1"
    brucella_identifications["1111101111111111"] = "suis2"
    brucella_identifications["1111110111111101"] = "suis3"
    brucella_identifications["1111111011111111"] = "ceti1"
    brucella_identifications["1111111001111111"] = "ceti1"
    brucella_identifications["1111111101111111"] = "ceti2"
    brucella_identifications["1111111110111101"] = "suis4"
    brucella_identifications["1111111110011101"] = "canis"
    brucella_identifications["1111111111101111"] = "ovis"

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
    bovis_identifications["11001110"] = "af" #bovis
    bovis_identifications["11011110"] = "af" #bovis
    bovis_identifications["11001100"] = "af" #bovis
    
    para_identifications = {}
    para_identifications["1"] = "para"
    para_identifications["01"] = "para"
    para_identifications["11"] = "para"

    count_summary={}

#    with futures.ProcessPoolExecutor() as pool:
#        for v, count in pool.map(finding_best_ref, oligo_dictionary.values()):
#            for k, value in oligo_dictionary.items():
#                if v == value:
#                    count_summary.update({k:count})
#                    count_summary=OrderedDict(sorted(count_summary.items()))

    pool = multiprocessing.Pool()
    func = partial(finding_best_ref, R1unzip, R2unzip)
    for v, count in pool.map(func, oligo_dictionary.values()):
        for k, value in oligo_dictionary.items():
            if v == value:
                count_summary.update({k:count})
                count_summary=OrderedDict(sorted(count_summary.items()))
    pool.close()
    pool.join()

    count_list=[]
    for v in count_summary.values():
        count_list.append(v)
    brucella_sum=sum(count_list[:16])
    bovis_sum=sum(count_list[16:24])
    para_sum=sum(count_list[24:])
    
    print("Best reference Brucella counts:", file=write_out)
    for i in count_list[:16]:
        print(i,  end=',', file=write_out)
    
    print("\nBest reference TB counts:", file=write_out)
    for i in count_list[16:24]:
        print(i,  end=',', file=write_out)

    print("\nBest reference Para counts:", file=write_out)
    for i in count_list[24:]:
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
    brucella_binary=binary_list[:16]
    brucella_string=''.join(str(e) for e in brucella_binary)
    bovis_binary=binary_list[16:24]
    bovis_string=''.join(str(e) for e in bovis_binary)
    para_binary=binary_list[24:]
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










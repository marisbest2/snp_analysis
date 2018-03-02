from Bio.SeqIO.QualityIO import FastqGeneralIterator
from concurrent import futures
import regex
from collections import OrderedDict

def finding_sp(v):
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

def spoligo(R1un, R2un):

    global R1unzip
    global R2unzip
    R1unzip = R1un
    R2unzip = R2un
    
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
    
#max_workers=limited_cpu_count

    with futures.ProcessPoolExecutor() as pool:
        for v, count in pool.map(finding_sp, spoligo_dictionary.values()):
            for k, value in spoligo_dictionary.items():
                if v == value:
                    count_summary.update({k:count})
                    count_summary=OrderedDict(sorted(count_summary.items()))
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

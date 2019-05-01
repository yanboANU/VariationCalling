import os
import sys
#from itertools import ifilter,imap
import collections
import string


MAXINT = 600000000

######################################
#read one list, only contain one list
#one line by one line
######################################
def read_string_list(filename): # read_file_list
    f = open(filename, "r")
    fileList = []
    for line in f:
        fileList.append(line.strip())
    return fileList 


def read_int_list(filename): # read_mutation_list
    f = open(filename, "r")
    pos = []  
    for line in f:  
        words = line.split()
        #if len(words) == 1 or len(words) == 3 or len(words) == 5:
        if words[0].isdigit():
            pos.append(int(words[0]))
    return pos

#this list must in order
def read_int_list_local(filename, start, end, base): #read_multiple_pos
    multiplePos = set() 
    f = open(filename, "r")
    for line in f:
        words = line.strip().split(" ")
        val = int(words[0])  
        if val < start:
            continue
        if val > end:
            break
        multiplePos.add(val+base)
    return multiplePos   



######################################
#read two lists, format a map, key and value both int
#one line by one line
######################################
def read_cov(filename):
    coverage = {}
    f = open(filename, "r")
    count = 0
    for line in f:
        if line.startswith("p"):
            continue
        words = line.strip().split(" ")
        coverage[int(words[0])] = int(words[1])
        if int(words[1]) >= 8:
            count += 1
    return coverage, count



#####################################
#read snp, delete 
#
#####################################

def read_pre_snp(filename):  #def read_snp(filename):
    f = open(filename, "r")
    snpPosition = set()
    snpContent = {}
    snpSupportNum = {}
    for line in f: 
        #print line
        if line.startswith('#'):
            continue
        words = line.strip().split()
        if len(words) == 5: #for file contig_snp_mutation
            pos = int(words[0])
            snpPosition.add(pos) 
            snpContent[pos] = (words[1], words[3])
            snpSupportNum[pos] = (int(words[2]), int(words[4]))
    return snpPosition, snpContent, snpSupportNum

def read_real_snp(filename, start=0, end=MAXINT, base=0):
# read_snp2
    f = open(filename, "r")
    snpPosition = set()
    snpContent = {}
    #requre snp store in order
    count =0
    for line in f:  
        words = line.strip().split()
        count += 1
        if len(words) == 4 or len(words) == 3: #for file mutation record
            pos = int(words[0])
            if pos < start:
                continue
            if pos > end:
                break
            snpPosition.add(pos+base)
            snpContent[pos+base] = (words[1], words[2])
        else:
            print ("error", line)
    print (count)       
    f.close()
    return snpPosition, snpContent

def read_delete(filename, start=0, end=MAXINT, base=0):

    f = open(filename, "r")
    snpPosition = set()
    snpContent = {}
    #requre snp store in order
    for line in f:  
        words = line.strip().split()
        if len(words) == 4 or len(words) == 3: #for file mutation record
            pos = int(words[0])
            if pos < start:
                continue
            if pos > end:
                break
            snpPosition.add(pos+base)
            snpContent[pos+base] = (words[1], words[2])
            deleteLen = len(words[1])-1
            count = 1
            if deleteLen > 1:
                snpPosition.add(pos+base+count)
                count += 1
                deleteLen -= 1
    return snpPosition, snpContent








###################################
#real two line, all number
#ref: 100, 101, 102, 103
#read:  0,   1,   2,   3
###################################
def read_align(filename):
    f = open(filename, "r")
    ref_pos = []
    contig_pos = []
    ref_label = 1 
    for line in f:
        words = line.split(',')
        if words[0].isdigit():
            if len(ref_pos) == 0:
                for c in words:
                    ref_pos.append(int(c))   
            else: 
                for c in words:
                    contig_pos.append(int(c)) 
    assert len(ref_pos) == len(contig_pos)
    '''
    fout = open("ref_contig", "w")  
    for i in range(len(ref_pos)):
        fout.write("%s %s\n" % (ref_pos[i], contig_pos[i]))
    '''
    return ref_pos, contig_pos



def read_phasing_result(filename):
    ##########
    #unfinish only read one segment
    ##########
    f = open(filename, "r")
    haplos = [] 
    lineNumber = 0
    for line in f:
        haplotype = {}
        lineNumber += 1
        if lineNumber % 6 == 0: 
            words = line.strip().split(',')
            #print (words)
            #print (len(words))    
            continue
        elif lineNumber % 6 == 1 and lineNumber >1:
            binarySeq = line.strip()
        else:
            continue 
        print ("len binnary",len(binarySeq))
        print ("len position",len(words))
        assert len(binarySeq) == len(words)
        for i in range(len(words)):
            haplotype[int(words[i])] = binarySeq[i] 
        haplos.append(haplotype) 
    return haplos

def read_blasr_m5(fileName):
    f = open(fileName,"r")
    for line in ifilter(lambda x: len(x) > 0, imap(string.strip, f)):
        words = line.split()
        if len(words) <= 4:
            continue
        #print words[:4] 
        queryName, queryLen, queryS, queryE = words[:4]
        queryDirection = words[4]
        targetName, targetLen, targetS, targetE = words[5:5+4]
        targetDirection = words[9] 
 
        querySeq = words[-3]   # ATCG-
        align = words[-2]      # |,*
        targetSeq = words[-1]  # ATCG-    
 
    f.close()
    query = sequence.Sequence(queryName, queryLen, queryS, queryE, querySeq, queryDirection) 
    target = sequence.Sequence(targetName, targetLen, targetS, targetE, targetSeq, targetDirection)
    alignObj = alignment.Alignment(query, target, align)  # blasr first reference, second contig
    #alignObj = alignment.Alignment(target, query, align)  # blasr first contig, second reference
    return alignObj



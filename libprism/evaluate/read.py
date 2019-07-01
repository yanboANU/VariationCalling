#########################################################################
# File Name: read.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Thu 30 May 2019 15:34:53 AEST
#########################################################################
#!/bin/bash


def read_vcf(filename):
    f = open(filename, "r")
    snps = {}
    ID = 1
    for line in f:
        if line.startswith('#'):
            continue
        words = line.split()
        pos = int(words[1])
        #print (pos)
        s1 = words[3].strip()
        s2 = words[4].strip()
        homo = words[9].split(':')[0]
        #print homo
        assert homo == "1/0" or homo == "1|0" or homo == "0|1" or homo == "0/1"
        s1Len = len(s1)
        s2Len = len(s2)
        assert s1Len == 1 and s2Len == 1
        if s1 == '.' or s2 == '.':
            continue
        snps[pos] = (s1, s2, ID)
        ID += 1

    f.close()     
    return snps    


def read_group_kmer(filename):
    pair_kmer = []
    pair_group_length = []
    state = 0
    with open(filename, "r") as f:
        for line in f:
            words= line.strip().split()
            if line.startswith('#'):
                continue
            if state == 0 and line.startswith("group"):
                pair_kmer.append(words[1:])
                state = 1
            elif state == 1:
                a = len(words)
                state = 2
            elif state == 2:
                b = len(words)
                state = 0
                pair_group_length.append( (a, b ) )

    return pair_kmer, pair_group_length      

def read_multip_columns(filename):
    pair_kmer = []
    with open(filename, "r") as f:
        for line in f:
            words= line.strip().split()
            pair_kmer.append(words)
    return pair_kmer      


def read2columns(filename):
    ID2ID = {}
    with open(filename, "r") as f:
        for line in f:
            words = line.strip().split()
            ID2ID[ int(words[0]) ] = int(words[1])
    
    return ID2ID       

def read1column(filename, pos):

    #ID = set()
    ID = list()
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            words = line.split()
            #ID.add( int( words[pos] ) )
            ID.append( int( words[pos] ))
    return ID

def get_diff(filename, pos):
    nums = []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            words = line.split()
            nums.append( int( words[pos] ) )
    numsLen = len(nums)
    diff = []
    for i in range(numsLen-1):
        diff.append(nums[i+1]-nums[i])
    diffSorted = sorted(diff)    
    for ele in diffSorted:
        print (ele)



def read1column_str(filename):

    ID = set()
    with open(filename, "r") as f:
        for line in f:
            words = line.split()
            ID.add( words[0] )
    return ID


def read_matrix(filename):
    groupID = {}
    readID = 1
    reads = {}
    with open(filename, "r") as f:
        for line in f:
            if readID % 10000 == 0:
                print ("deal read number", readID)
            words = line.split()
            num = int(words[0])
            reads [ readID ] = {}
            for i in range(num):
                ID = int(words[2*i+2])
                reads[readID][ID] = int(words[2*i+3])
                if ID not in groupID:
                    groupID[ID] = set()
                groupID[ID].add(readID)
            readID += 1
    return groupID, reads     
            
#get_diff("SNP_pos", 0)

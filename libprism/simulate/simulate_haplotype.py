#########################################################################
# File Name: get_2_haplotype.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Mon 24 Jun 2019 14:15:36 AEST
#########################################################################
#!/bin/bash

import collections
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import os
import sys

##########################
#input: ref and vcf 
#output: two haplotypes
##########################
#vcfFile = "/home/yulin/software/HapCUT2/reproduce_hapcut2_paper/run_hapcut2_fosmid/data/NA12878_hg19_VCFs/chr15.hg19.vcf"
snpFile = "/home/yulin/bio/VariationCalling/data/NA12878/VCF/NA12878_hg18_snp_VCFs/" + sys.argv[1] + ".vcf" 
indelFile = "/home/yulin/bio/VariationCalling/data/NA12878/VCF/NA12878_hg18_indel_VCFs/" + sys.argv[1] + ".vcf"
refFile = "/home/yulin/bio/Data/reference/NCBI36_hg18/" + sys.argv[1] + ".fa"
outFile = "/home/yulin/bio/VariationCalling/data/NA12878/reference/" + sys.argv[1] + "_2_haplotypes.fa"

def read_vcf(filename, vcfType):
    # vcfType=0 -> SNP; vcfType=1 -> indel
    mutations = []
    with open(filename, "r")as f:
        for line in f:
            if line.startswith('#'):
                continue
            words = line.split()
            pos = int(words[1])
            s1 = words[3].strip()
            s2 = words[4].strip()
            homo = words[9].split(':')[0]
            #print homo
            assert homo == "1/0" or homo == "1|0" or homo == "0|1" or homo == "0/1"
            s1Len = len(s1)
            s2Len = len(s2)
            if vcfType == 0:
                assert s1Len == 1 and s2Len == 1
            elif vcfType == 1:
                assert s1Len > 1 or s2Len > 1
            if s1 == '.' or s2 == '.':
                continue
            mutations.append( (pos, s1, s2, homo) )
    return mutations    


def simulate_haplotypes_only_snp(refFile, snpFile):

    snps = read_vcf(snpFile, 0)
    print ("finish read snps")
    record = SeqIO.read(open(refFile), "fasta")
    seqH1List = list(record.seq)
    seqH2List = list(record.seq)
    for (pos, s1, s2, homo) in snps:
        #print pos, s1, s2
        assert record.seq[pos-1] == s1 or record.seq[pos-1].upper() == s1 
        if homo == "1|0" or homo == "1/0":
            seqH1List[pos-1] = s1
            seqH2List[pos-1] = s2
        elif homo == "0|1" or homo == "0/1":
            seqH1List[pos-1] = s2
            seqH2List[pos-1] = s1

    seqH1 = ''.join(ele for ele in seqH1List)
    seqH2 = ''.join(ele for ele in seqH2List)
    records = []
    rec1 = SeqRecord( Seq(seqH1) , id="chr15_H1")
    rec2 = SeqRecord( Seq(seqH2) , id="chr15_H2")
    records.append(rec1)
    records.append(rec2)
    SeqIO.write(records, outFile, 'fasta') 

    # for test
    '''
    test = snps[:10]
    for record in SeqIO.parse(outFile, "fasta"):
        print(record.id)       
        for (pos, s1, s2, homo) in test:
            print pos, record.seq[pos-1], s1, s2,
    '''
#chr1   35931175    .   AT  A
#chr1   35931176    rs11264189  T   A
#both in ground-truth,  
#def simulate_haplotypes_(refFile, snpFile, indelFile):
def simulate_haplotypes(refFile, snpFile, indelFile, indelLen=1): 
    # indelLen=-1, consider all indel in vcf, indelLen = 1, only consider indelLen = 1
    mutations = []
    snps = read_vcf(snpFile, 0)
    print ("finish read snps, snp number: %s" % len(snps))
    indels = read_vcf(indelFile, 1)
    print ("finish read indels, indel number: %s" % len(indels))
    mutations.extend(snps)
    mutations.extend(indels)
    sortedMutations = sorted(mutations) 
    record = SeqIO.read(open(refFile), "fasta")
    seqH1List = list(record.seq)
    seqH2List = list(record.seq)
    adjustPos1 = 0
    adjustPos2 = 0
    assert indelLen == 1
    snpCount, insertCount, delCount = 0, 0, 0
    prePos = -1
    for (pos, s1, s2, homo) in sortedMutations:
        #print pos, s1, s2
        if pos == prePos:
            continue
        assert ( record.seq[pos-1].upper() == s1[0:1] )
        if len(s1) == 1 and len(s2) == 1:
            snpCount += 1
            if homo == "1|0" or homo == "1/0":
                print adjustPos1, pos
                if seqH1List[pos + adjustPos1 -1].upper() != s1:
                    print prePos, pos, "beacuse of neigbor, do not add this mutation"
                    snpCount -= 1
                    continue
                assert(seqH1List[pos + adjustPos1 -1].upper() == s1)
                seqH2List[pos + adjustPos2 -1] = s2
            elif homo == "0|1" or homo == "0/1":
                assert(seqH2List[pos + adjustPos2 -1].upper() == s1)
                seqH1List[pos + adjustPos1 -1] = s2
            prePos = pos    
        elif len(s1) == 1 and len(s2) == 2: # insert
            insertCount += 1
            assert s1 == s2[0:1]
            if homo == "1|0" or homo == "1/0":
                assert(seqH1List[pos + adjustPos1 -1].upper() == s1)
                seqH2List.insert(pos + adjustPos2, s2[1:2])
                adjustPos2 += 1
                print "insert at H2"
                print "pos, adjustPos1, adjustPos2", pos, adjustPos1, adjustPos2
            elif homo == "0|1" or homo == "0/1":
                assert(seqH2List[pos + adjustPos2 -1].upper() == s1)
                seqH1List.insert(pos + adjustPos1, s2[1:2] )
                adjustPos1 += 1
                print "insert at H1"
                print "pos, adjustPos1, adjustPos2", pos, adjustPos1, adjustPos2
            prePos = pos    
        elif len(s1) == 2 and len(s2) == 1: # delete
            delCount += 1
            assert s1[0:1] == s2
            if homo == "1|0" or homo == "1/0":
                assert( ''.join(ele.upper() for ele in seqH1List[pos + adjustPos1 -1 : pos+adjustPos1+1]) == s1 )
                seqH2List = seqH2List[:pos+adjustPos2 ] + seqH2List[pos+adjustPos2+1:]
                adjustPos2 -= 1
                print "delete at H2"
                print "pos, adjustPos1, adjustPos2", pos, adjustPos1, adjustPos2
            elif homo == "0|1" or homo == "0/1":
                #print seqH2List[pos + adjustPos2 -1 : pos+adjustPos2+1], s1
                #print ''.join(ele.upper() for ele in seqH2List[pos + adjustPos2 -1 : pos+adjustPos2+1])
                assert( ''.join(ele.upper() for ele in seqH2List[pos + adjustPos2 -1 : pos+adjustPos2+1]) == s1)
                seqH1List = seqH1List[:pos+adjustPos1 ] + seqH1List[pos+adjustPos1+1:]
                adjustPos1 -= 1
                print "delete at H1"
                print "pos, adjustPos1, adjustPos2", pos, adjustPos1, adjustPos2
            prePos = pos    

    print ("snp number: %s; insert number: %s; delete number: %s" % (snpCount, insertCount, delCount) )
    seqH1 = ''.join(ele for ele in seqH1List)
    seqH2 = ''.join(ele for ele in seqH2List)
    records = []
    rec1 = SeqRecord( Seq(seqH1) , id=sys.argv[1] + "_H1")
    rec2 = SeqRecord( Seq(seqH2) , id=sys.argv[1] + "_H2")
    records.append(rec1)
    records.append(rec2)
    SeqIO.write(records, outFile, 'fasta')

    
simulate_haplotypes(refFile, snpFile, indelFile) 

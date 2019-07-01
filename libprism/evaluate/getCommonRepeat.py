#########################################################################
# File Name: getCommonRepeat.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Tue 18 Jun 2019 09:49:33 AEST
#########################################################################
#!/bin/bash

import collections
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import sys
import tools

import read


lowNucleotide = set('atcg')

refFilename ="/media/yanbo/Data/reference/hg37/chr" + sys.argv[1] + ".fa"

record = SeqIO.read(open(refFilename), "fasta")
print (record.id)
seq = str(record.seq)
seqLen = len(seq)

print ("chr", sys.argv[1], "length:", seqLen)
print ("none N length:", seqLen-seq.count('N') )
commonRepeatRegion = []
for i in range(seqLen):
    if seq[i] in lowNucleotide:
        commonRepeatRegion.append(i)

print ("common repeat region length", len(commonRepeatRegion) )
for ele in commonRepeatRegion:
    print (ele)

vcfFilename ="/home/yanbo/software/HapCUT2/reproduce_hapcut2_paper/run_hapcut2_fosmid/data/NA12878_hg19_VCFs/chr" + sys.argv[1] + ".hg19.vcf"
snps = read.read_vcf(vcfFilename)
snpsList = sorted( list(snps.keys()) )


lenSNP, lenR = len(snpsList), len(commonRepeatRegion)


#fout = open("GroupID2RefPos2", "w")
inter = []
indexSNP, indexR = 0, 0
count = 0 
while indexSNP < lenSNP and indexR < lenR:
    if (snpsList[indexSNP] < commonRepeatRegion[indexR]):
        indexSNP +=1
    elif snpsList[indexSNP] > commonRepeatRegion[indexR]:   
        indexR +=1
    elif snpsList[indexSNP] == commonRepeatRegion[indexR]:
        temp = snpsList[indexSNP]
        tempSeq = seq[temp-25:temp+26]
        #print (temp, tempSeq)
        if tempSeq.count('A') + tempSeq.count('T') + tempSeq.count('C') + tempSeq.count('G') == 0:
            count += 1
        inter.append( commonRepeatRegion[indexR] )
        indexSNP +=1
        indexR +=1


print ("true SNP in common repeat region", len(inter) )
print ("all 51mer of SNP in common repeat region", count )

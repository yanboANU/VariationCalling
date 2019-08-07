#########################################################################
# File Name: transfer_discoSNP_output.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Wed Jul  3 19:44:14 2019
#########################################################################
#!/bin/bash

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import os
import sys
import tools
#######################################
#input: .snp 
#output: snp_pair_kmer indel_pair_kmer
#######################################
#vcfFile = "/home/yulin/software/HapCUT2/reproduce_hapcut2_paper/run_hapcut2_fosmid/data/NA12878_hg19_VCFs/chr15.hg19.vcf"

inFile = sys.argv[1]
filterCov = sys.argv[2]
#inFile = "output.5.snp"

snpFile = "snp_31mer_cov" + filterCov + "_pair_kmer"
indelFile = "indel_31mer_cov"+ filterCov +"_pair_kmer"

print ("input: %s" % (inFile) )
# python3 sorted snp pair kmer correct
# python2 sorted snp pair kmer fail

snpPairKmer, indelPairKmer = [], []
state = 0
snpNumber = 0
with open(inFile, "r") as f:
    line = f.readline()
    while line:
        if state == 0 and line.startswith(">"):
            words = line.split("_")
            ID = words[0][1:]
            mutation = words[5]
            if mutation == "INDEL":
                indelLen = len(line.strip().split(":")[-1]) - 1
            state = 1 
            line = f.readline()
            continue
        elif state == 1:
            kmer1 = line.strip()
            state = 2
            line = f.readline()
            line = f.readline()
            continue
        elif state == 2:
            kmer2 = line.strip()
            if mutation == "SNP":
                snpNumber += 1 
                #print (ID, kmer1, kmer2)
                assert len(kmer1) == len(kmer2)
                #kmer1 = kmer1[15:-15] # if we want middle 31mer
                #kmer2 = kmer2[15:-15] # this two line 
                if tools.hamming_distance(kmer1, kmer2) != 1:
                    #print ("not isolated with 31mer")
                    #print (ID, kmer1, kmer2)
                    state = 0
                    line = f.readline()
                    continue
                smallerKmer1, smallerKmer2 = tools.get_smaller_pair_kmer(kmer1, kmer2)
                snpPairKmer.append( (smallerKmer1, smallerKmer2) )
                #print (ID, smallerKmer1, smallerKmer2)
            elif mutation == "INDEL" and indelLen == 1:
                kmer2 = kmer2[1:]
                
                #kmer1 = kmer1[15:-15] # for middle 31mer
                #kmer2 = kmer2[15:-15]
                smallerKmer1, smallerKmer2 = tools.get_smaller_pair_kmer(kmer1, kmer2)
                indelPairKmer.append( (smallerKmer1, smallerKmer2) )
            state = 0
            line = f.readline()
        else:
            print ("error")

#print snpPairKmer
#print indelPairKmer

print ( "total snp number:", snpNumber )
print ( "isolated snp number:", len(snpPairKmer) )
sortedSnpPairKmer = sorted(snpPairKmer)
snpOUT = open(snpFile, "w")

#for (kmer1, kmer2, ID, pos) in sortedSnpPairKmer:
    #snpOUT.write("%s %s %s %s\n" % (kmer1, kmer2, ID, pos) )

for (kmer1, kmer2) in sortedSnpPairKmer:
    snpOUT.write("%s %s\n" % (kmer1, kmer2) )
snpOUT.close()


print ( "indel number", len(indelPairKmer) )
sortedIndelPairKmer = sorted(indelPairKmer)
indelOUT = open(indelFile, "w")
for (kmer1, kmer2) in sortedIndelPairKmer:
    indelOUT.write("%s %s\n" % (kmer1, kmer2) )
indelOUT.close()




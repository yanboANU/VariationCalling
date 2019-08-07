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
#input: .vcf and .fa
#output: snp_pair_kmer indel_pair_kmer
#######################################
#vcfFile = "/home/yulin/software/HapCUT2/reproduce_hapcut2_paper/run_hapcut2_fosmid/data/NA12878_hg19_VCFs/chr15.hg19.vcf"
vcfFile = sys.argv[1]
faFile = sys.argv[2] 

#vcfFile = "discoRes_k_31_c_3_D_100_P_3_b_0_coherent.vcf"
#faFile = "discoRes_k_31_c_3_D_100_P_3_b_0_coherent.fa"

snpFile = "snp_pair_kmer"
indelFile = "indel_pair_kmer"

print ("input: %s %s" % (vcfFile, faFile) )
# python3 sorted snp pair kmer correct
# python2 sorted snp pair kmer fail

def read_vcf(filename):
    mutations = {}
    with open(filename, "r")as f:
        for line in f:
            if line.startswith('#'):
                continue
            words = line.split()
            contigID = words[0] 
            pos = int(words[1])
            s1 = words[3].strip()
            s2 = words[4].strip()
            homo = words[9].split(':')[0]
            #print homo
            if homo == "0/0" or homo == "0|0" or homo == "1|1" or homo == "1/1":
                continue
            assert homo == "1/0" or homo == "1|0" or homo == "0|1" or homo == "0/1"
            s1Len = len(s1)
            s2Len = len(s2)
            if contigID not in mutations:
                mutations[contigID] =[]
            mutations[contigID].append( ( pos, s1, s2, homo) )
    return mutations    

mutations = read_vcf(vcfFile)
print ( "mutations number:", len(mutations) )
L =15

snpPairKmer = []
indelPairKmer = []
for record in SeqIO.parse(faFile, "fasta"):
    ID = record.id.split('|')[0] 
    if ID in mutations:
        seq = record.seq.upper()
        mLen = len(mutations[ID])
        for (pos, s1, s2, homo) in mutations[ID]:
            if len(s1) + len(s2) >3:
                continue
            if mLen > 1:
                pos = pos - 1
                # one mutation, vcf position starts from 0
                # more than on mutation, vcf position starts from 1  
            if len(s1) == 1 and len(s2) == 1:
                if seq[pos] != s1.upper():
                    print ("snp error pos range", ID, pos, seq[pos-2:pos+3], s1, s2)
                    print (seq)
                    continue
                assert seq[pos] == s1.upper()
                kmer1 = seq[pos-15:pos+16]
                kmer2 = seq[pos-15:pos] + s2.upper() + seq[pos+1:pos+16]
                assert len(kmer1) + len(kmer2) == 62
                smallerKmer1, smallerKmer2 = tools.get_smaller_pair_kmer(kmer1, kmer2)
                #snpPairKmer.append( (smallerKmer1, smallerKmer2, ID, pos) )        
                snpPairKmer.append( (smallerKmer1, smallerKmer2) )    
            elif len(s1) == 1 and len(s2) == 2:
                #print "indel success", ID, pos, seq[pos-2:pos+3], s1, s2 
                #print seq
                print ("not happen") #  assert s1[0:1] == s2[0:1]
                kmer1 = seq[pos-14:pos+16]
                kmer2 = seq[pos-14:pos] + s2.upper() + seq[pos+1:pos+16]
                assert len(kmer1) + len(kmer2) == 61
                smallerKmer1, smallerKmer2 = tools.get_smaller_pair_kmer(kmer1, kmer2)
                indelPairKmer.append( (smallerKmer1, smallerKmer2) )
            elif len(s1) == 2 and len(s2) == 1:   # vcf pos record correct but s1 and s2 not correct
                #print "indel success", ID, pos, seq[pos-2:pos+3], s1, s2 
                #print seq  #assert s1[0:1] == s2[0:1]
                kmer1 = seq[pos-14:pos+17]
                kmer2 = seq[pos-14:pos+1] + seq[pos+2:pos+17]
                assert len(kmer1) + len(kmer2) == 61
                smallerKmer1, smallerKmer2 = tools.get_smaller_pair_kmer(kmer1, kmer2)
                indelPairKmer.append( (smallerKmer1, smallerKmer2) )

print ( "snp number:", len(snpPairKmer) )
sortedSnpPairKmer = sorted(snpPairKmer)
print ( "sorted snp number", len(sortedSnpPairKmer) )
snpOUT = open(snpFile, "w")

#for (kmer1, kmer2, ID, pos) in sortedSnpPairKmer:
    #snpOUT.write("%s %s %s %s\n" % (kmer1, kmer2, ID, pos) )

for (kmer1, kmer2) in sortedSnpPairKmer:
    snpOUT.write("%s %s\n" % (kmer1, kmer2) )
snpOUT.close()


sortedIndelPairKmer = sorted(indelPairKmer)
indelOUT = open(indelFile, "w")
for (kmer1, kmer2) in sortedIndelPairKmer:
    indelOUT.write("%s %s\n" % (kmer1, kmer2) )
indelOUT.close()




#########################################################################
# File Name: getKmerFromVCF_REF.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Thu 09 May 2019 10:45:06 AEST
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




def write_pair_kmer(outFile, kmers):

    sortedKmers = sorted(kmers)
    with open(outFile, "w") as f:
        for (kmer1, kmer2) in sortedKmers:
            #f.write("%s %s %s %s %s %s %s %s\n" % (ele[0], ele[1], ele[2], ele[3], ele[4], ele[5], tools.reverse(ele[0]), tools.reverse(ele[1]) ) )
            f.write("%s %s\n" % (kmer1, kmer2) )


def get_snp_pair_kmer(vcfFilename):

    snps = read.read_vcf(vcfFilename)
    kmerFilename="chr" + sys.argv[1] + ".snp.real." + sys.argv[2] + "mer"
    kmers = []
    for key in snps:
        
        #assert seq[key-1] == snps[key][0] or seq[key-1] == snps[key][1]
        assert seq[key-1] == snps[key][0] 
        h1 = seq[key-int(k/2)-1 : key-1] + snps[key][0] + seq[key : key+int(k/2) ] # 0
        h2 = seq[key-int(k/2)-1 : key-1] + snps[key][1] + seq[key : key+int(k/2) ] # 1
        if h1.count('N') > 0 or h2.count('N') > 0:
            continue
        '''
        new_h1 = tools.reverse(h1) # 0
        new_h2 = tools.reverse(h2) # 1
        min_h= min(h1,h2)
        min_newh = min(new_h1, new_h2)
        ID = snps[key][2]
        if min_h < min_newh: 
            if h1 < h2:
                kmers.append( (h1, h2, key, ID, 0, 1) )
            else:
                kmers.append( (h2, h1, key, ID, 1, 0) )
        else:        
            if new_h1 < new_h2:
                kmers.append( (new_h1, new_h2, key, ID, 0, 1) )
            else:
                kmers.append( (new_h2, new_h1, key, ID, 1, 0) )
        '''
        smallerH1, smallerH2 = tools.get_smaller_pair_kmer(h1, h2)
        kmers.append( (smallerH1, smallerH2) )

    write_pair_kmer(kmerFilename, kmers)


def get_indel_pair_kmer(vcfFilename):

    indels = read.read_vcf(vcfFilename)
    kmerFilename="chr" + sys.argv[1] + ".indel.real." + sys.argv[2] + "mer"
    kmers = []
    for key in indels:
        s1, s2, ID = indels[key]
        lenS1, lenS2 = len(s1), len(s2)
        if lenS1 + lenS2 > 3:
            continue
        assert lenS1 + lenS2 >= 2
        if len(s1) == 1 and len(s2) == 2:
            assert seq[key-1] == s1
            #print s1, s2, seq[key-2:key]
            assert s2[0] != s2[1]
            #while s2[1] == seq[key-1]: # delete content is s2[1]
                #key = key-1            # delete happen at "AAA" region, always think delete first poisition                  
            h1 = seq[key-int(k/2) : key+int(k/2)] # k-1 
            h2 = seq[key-int(k/2) : key-1] + s2 + seq[key : key+int(k/2)] # 1 # len: k
            #print "aa", len(h1), len(h2)
            assert len(h1) == k-1 and len(h2) == k
            h1, h2 = h2, h1

            if h1.count('N') > 0 or h2.count('N') > 0:
                continue
        elif len(s1) == 2 and len(s2) == 1:   
            assert seq[key-1:key+1] == s1
            assert s1[0] != s1[1]            
            #while s1[1] == seq[key-1]:
                #key = key-1
            h1 = seq[key-int(k/2) : key+int(k/2)+1] # k 
            h2 = seq[key-int(k/2) : key] + seq[key+1 : key+int(k/2)+1] # 1 # len: k-1
            #print "bb", len(h1), len(h2)
            assert len(h1) == k and len(h2) == k-1
            if h1.count('N') > 0 or h2.count('N') > 0:
                continue
        smallerH1, smallerH2 = tools.get_smaller_pair_kmer(h1, h2)
        kmers.append( (smallerH1, smallerH2) )

    write_pair_kmer(kmerFilename, kmers)


'''
allFile = "chr" + sys.argv[1] + ".all." + sys.argv[2] + "mer"
foutAll = open(allFile, "w")
for i in range(seqLen-21):
    mer = seq[i:i+k]
    if mer.count('N') > 0:
        continue
    Rmer = tools.reverse(mer)
    if Rmer < mer:
        mer = Rmer
    foutAll.write("%s %s\n" % (mer, i))    
foutAll.close()       
'''


# this simulate data is based on hg18
refFilename="/home/yulin/bio/Data/reference/NCBI36_hg18/chr22.fa"
snpVCFFile="/home/yulin/bio/VariationCalling/data/NA12878/VCF/NA12878_hg18_snp_VCFs/chr22.vcf"
indelVCFFile="/home/yulin/bio/VariationCalling/data/NA12878/VCF/NA12878_hg18_indel_VCFs/chr22.vcf"
        
# this illumina data align to hg19
#print ("input chrID kmer-size")
#refFilename ="/home/yulin/bio/Data/reference/GRCh37_hg19/chr" + sys.argv[1] + ".fa"
#vcfFilename ="/home/yulin/software/HapCUT2/reproduce_hapcut2_paper/run_hapcut2_fosmid/data/NA12878_hg19_VCFs/chr" + sys.argv[1] + ".phased.vcf"


record = SeqIO.read(open(refFilename), "fasta")
print (record.id)
seq = str(record.seq).upper()
seqLen = len(seq)
k=int(sys.argv[2])
get_snp_pair_kmer(snpVCFFile)
get_indel_pair_kmer(indelVCFFile)


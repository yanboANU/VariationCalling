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
        
# this illumina data align to hg19
print ("input chrID kmer-size")
refFilename ="/media/yanbo/Data/reference/hg37/chr" + sys.argv[1] + ".fa"
vcfFilename ="/home/yanbo/software/HapCUT2/reproduce_hapcut2_paper/run_hapcut2_fosmid/data/NA12878_hg19_VCFs/chr" + sys.argv[1] + ".hg19.vcf"

#record = SeqIO.read(open(sys.argv[1]), "fasta")
record = SeqIO.read(open(refFilename), "fasta")
print (record.id)
seq = str(record.seq).upper()
seqLen = len(seq)
snps = read.read_vcf(vcfFilename)
k=int(sys.argv[2])

kmerFilename="chr" + sys.argv[1] + ".real." + sys.argv[2] + "mer"
kmers = []

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


for key in snps:
    
    #print (seq[key-1], snps[key][0], snps[key][1])
    assert seq[key-1] == snps[key][0] or seq[key-1] == snps[key][1]
    h1 = seq[key-int(k/2)-1 : key-1] + snps[key][0] + seq[key : key+int(k/2) ] # 0
    h2 = seq[key-int(k/2)-1 : key-1] + snps[key][1] + seq[key : key+int(k/2) ] # 1
    if h1.count('N') > 0 or h2.count('N') > 0:
        continue
    #print (h1)
    #print (h2)
    #print (seq[key-int(k/2)-1 : key+int(k/2) ] , key)

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

sortedKmers = sorted(kmers)
with open(kmerFilename, "w") as f:
    for ele in sortedKmers:
        f.write("%s %s %s %s %s %s %s %s\n" % (ele[0], ele[1], ele[2], ele[3], ele[4], ele[5], tools.reverse(ele[0]), tools.reverse(ele[1]) ) )



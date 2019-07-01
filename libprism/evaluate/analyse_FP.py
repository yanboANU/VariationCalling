#########################################################################
# File Name: compareRealRair.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Thu 09 May 2019 11:51:38 AEST
#########################################################################
#!/bin/bash
import sys
import tools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

def read_pair_kmer(filename):
    pair_kmer = []
    with open(filename, "r") as f:
        for line in f:
            words = line.strip().split()
            pair_kmer.append(words)
    return pair_kmer      



def get_FP_position(k1, k2):
    refFilename ="/media/yanbo/Data/reference/hg37/chr22.fa"
    record = SeqIO.read(open(refFilename), "fasta")
    print (record.id)
    seq = str(record.seq).upper()
    #count = 0
    #for (k1, k2, c1, c2) in pairFP:
    Rk1=tools.reverse(k1)
    Rk2=tools.reverse(k2) 
    a,b = seq.count(k1), seq.count(Rk1) 
    if a >0:
        print (k1, k2, a, seq.index(k1) )
    if b >0:    
        print (Rk1, Rk2, b, seq.index(Rk1) )
    
    c,d = seq.count(k2), seq.count(Rk2) 
    if c >0:
        print (k2, k1, c, seq.index(k2) )
    if d>0:    
        print (Rk2, Rk1, d, seq.index(Rk2) )
    #if count >=10:
        #sys.exit()
    #count += 1



FPfile="FP"

pairFP = read_pair_kmer(FPfile)
count = 0
print (len(pairFP))
for (k1, k2, c1, c2) in pairFP:
    l = int(len(k1)/2)
    ele1=k1[0]
    ele2=k2[0]
    pre1 = k1[0:l]
    pre2 = k2[0:l]
    if pre1.count(ele1) == l  or pre2.count(ele2) == l:
        print ("pre", k1, k2)
        count += 1
    else: 
        ele1=k1[-1]
        ele2=k2[-1]
        sub1 = k1[-l:]
        sub2 = k2[-l:]
        if sub1.count(ele1) == l  or sub2.count(ele2) == l:
            print ("sub", k1, k2)
            count += 1
        else:
            print ("unknow", k1, k2, c1, c2)
            get_FP_position(k1, k2)
print (count)


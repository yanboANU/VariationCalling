#########################################################################
# File Name: variationCalling.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Mon 18 Mar 2019 13:59:49 AEDT
#########################################################################
#!/bin/bash
import argparse
import os
import sys
from libprismv2.local import kmercalling
#from libprism.local.prepare import clouds_from_refhap, merge_clouds, print_clouds
from libprism.evaluate import read
#from math import log, exp
#from pysam
import logging
import time

#############################################################################

parser = argparse.ArgumentParser(description="Variation calling for Single-individual based on NGS reads")

#parser.add_argument('--bam', help='path to alignment bam file', required=True)
#parser.add_argument('--ref', help='reference or scaffolds', required=True)
parser.add_argument('--t1', help='kmer freq txt file (dsk result)', required=True)
parser.add_argument('--t2', help='k-1mer freq txt file (dsk result)', required=True)
parser.add_argument('--c1', help='low coverage', required=True)
parser.add_argument('--c2', help='high coverage', required=True)
parser.add_argument('--k', help='kmer size', required=True)
parser.add_argument('--b', help='for simulate data 0, real data 1, can dicide indel threshold', required=True)
#parser.add_argument('--o', help='output file prefix', required=True)

args = parser.parse_args()
logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

#############################################################################
time1 = time.clock()
lowCov, highCov = int(args.c1), int(args.c2)

## try to filter more erronoure kmer
#if lowCov < 10:
#    lowCov = 10

# try this
lowCov = lowCov - 1
highCov = highCov + 1

k = int(args.k)

kmerCov = kmercalling.pick_smaller_unique_kmer( args.t1, 
             lowCov, highCov)
k_1merCov = kmercalling.pick_smaller_unique_kmer( args.t2, 
             lowCov, highCov)

time2 = time.clock()
# uniq.kmer
#uniqKmer = read.read_multip_columns(args.t1) 
#uniqK_1mer = read.read_multip_columns(args.t2)
logging.info( "uniq Kmer size: %s" % len(kmerCov) )
logging.info( "uniq (K-1)mer size: %s" % len(k_1merCov) )
logging.info( "finish reading (kmer cov) and (k-1mer cov) file, cost %s seconds" % (time2-time1))
left_index, right_index = kmercalling.build_left_right_kmer_index(kmerCov)

mapK, snpPair, leftKmer = kmercalling.find_snp_pair_kmer(kmerCov, k, left_index, right_index)
time3 = time.clock()
logging.info( "finish find snp, cost %s second" % (time3 - time2) )
logging.info( "finish find snp, leftKmer size %s" % len(leftKmer) )
nonPair = []
#nonPair, highRepeat = kmercalling.find_non_pair_kmer(uniqKmer, k, left_index, right_index)
nonPair, highRepeat, leftKmer = kmercalling.find_non_pair_kmer(leftKmer, kmerCov, k, left_index, right_index)
time4 = time.clock()
logging.info( "finish find non snp, cost %s second" % (time4 - time3) )
logging.info( "finish find non snp, leftKmer size %s" % len(leftKmer) )
#time4 = time.clock()
#adjustKmer = read.read_multip_columns("chr22_k31.adjust.sort.txt")
#highRepeat = set()
#left_index, right_index = kmercalling.build_left_right_kmer_index(adjustKmer)
indelPair = []
if args.b == "0":
    threshold = 0
else:
    threshold = int(k/4)
indelPair = kmercalling.find_indel_pair_kmer(leftKmer, kmerCov, k_1merCov, k, left_index, right_index, threshold)#, highRepeat)
time5 = time.clock()
logging.info( "finish find indel, cost %s seconds" % (time5 - time4) )
'''
hetePairs = []
f1 = "k_" + str(k)+"_pair.snp"
fout = open(f1, "w")

ID = 0
for (k1, k2, c1, c2) in snpPair:
    fout.write(">kmer_snp%s_1_cov_%s\n" % (ID, c1))
    fout.write("%s\n" % ( k1 ) )
    fout.write(">kmer_snp%s_2_cov_%s\n" % (ID, c2))
    fout.write("%s\n" % ( k2 ) )
    ID += 1
    hetePairs.add( (k1, k2, "snp"+str(ID) ) )

ID = 0
for (k1, k2, c1, c2) in indelPair:
    fout.write(">kmer_snp%s_1_cov_%s\n" % (ID, c1))
    fout.write("%s\n" % ( k1 ) )
    fout.write(">kmer_snp%s_2_cov_%s\n" % (ID, c2))
    fout.write("%s\n" % ( k2 ) )
    ID += 1
    hetePairs.add( (k1, k2, "indel"+str(ID) ) )


ID = 0
for eles in snpPair:
    kmers = []
    #print (eles)
    covs = "".join("_" + str(c) for (kmer,c) in eles)
    cnt = 1
    for (kmer, c) in eles:
        kmers.append(kmer)
        fout.write(">kmer_snp%s_%s_cov%s\n" % (ID, cnt, c))
        fout.write("%s\n" % ( kmer ) )
        cnt += 1
    ID += 1
    hetePairs.append( (kmers, "snp"+str(ID) + "cov" + covs) )

ID = 0
for eles in indelPair:
    kmers = []
    covs = "".join("_" + str(c) for (kmer,c) in eles)
    cnt = 1
    for (kmer, c) in eles:
        kmers.append(kmer)
        fout.write(">kmer_indels%s_%s_cov%s\n" % (ID, cnt, c))
        fout.write("%s\n" % ( kmer ) )
        cnt += 1
    ID += 1
    hetePairs.append( (kmers, "indels"+str(ID) + "cov" + covs) )


ID = 0
for (k1, k2, c1, c2, c3, c4) in nonPair:
    kmers = []
    fout.write(">kmer_nonsnp%s_1_cov_%s_cov%s\n" % (ID, c1, c3))
    fout.write("%s\n" % ( k1 ) )
    fout.write(">kmer_nonsnp%s_2_cov_%s_cov%s\n" % (ID, c2, c4))
    fout.write("%s\n" % ( k2 ) )
    ID += 1
    kmers.append(k1)
    kmers.append(k2)
    hetePairs.append( (kmers, "nonsnp"+str(ID) ) )

fout.close()

#print (type(k))

left_index, right_index = kmercalling.build_left_right_kmer_index(uniqKmer)
extendPair = kmercalling.extend_pair(hetePairs, left_index, right_index, k)


f2 = "k_" + str(k)+"_extend_pair.snp"
with open(f2, "w") as fout:
    for (eles, ID) in extendPair:
        cnt = 1
        for kmer in eles:
            fout.write(">kmer_%s_%s\n" % (ID, cnt))
            fout.write("%s\n" % ( kmer ) )
            cnt += 1
'''            

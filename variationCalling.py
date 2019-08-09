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
from libprism.local import kmercalling
#from libprism.local.prepare import clouds_from_refhap, merge_clouds, print_clouds
from libprism.evaluate import read
#from math import log, exp
#from pysam
import logging

#############################################################################

parser = argparse.ArgumentParser(description="Variation calling for Single-individual based on NGS reads")

#parser.add_argument('--bam', help='path to alignment bam file', required=True)
#parser.add_argument('--ref', help='reference or scaffolds', required=True)
parser.add_argument('--t1', help='kmer freq txt file (dsk result)', required=True)
parser.add_argument('--t2', help='k-1mer freq txt file (dsk result)', required=True)
parser.add_argument('--c1', help='low coverage', required=True)
parser.add_argument('--c2', help='high coverage', required=True)
parser.add_argument('--k', help='kmer size', required=True)
#parser.add_argument('--o', help='output file prefix', required=True)

args = parser.parse_args()
logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

#############################################################################

lowCov, highCov = int(args.c1), int(args.c2)
k = int(args.k)

#uniqKmer = kmercalling.pick_smaller_unique_kmer( args.t1, 
#             lowCov, highCov)

#uniqK_1mer = kmercalling.pick_smaller_unique_kmer( args.t2, 
#             lowCov, highCov)

uniqKmer = read.read_multip_columns(args.t1) 
uniqK_1mer = read.read_multip_columns(args.t2)

logging.info( "finish reading (kmer cov) and (k-1mer cov) file" )

mapK, snpPair = kmercalling.find_snp_pair_kmer(uniqKmer, k)

logging.info( "finish find snp" )
indelPair = kmercalling.find_indel_pair_kmer(mapK, uniqK_1mer, k)

logging.info( "finish find indel" )
nonPair = kmercalling.find_non_pair_kmer(uniqKmer, k)

logging.info( "finish find non snp" )
hetePairs = set()
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
for (k1, k2, c1, c2, c3, c4) in nonPair:
    fout.write(">kmer_nonsnp%s_1_cov_%s_cov%s\n" % (ID, c1, c3))
    fout.write("%s\n" % ( k1 ) )
    fout.write(">kmer_nonsnp%s_2_cov_%s_cov%s\n" % (ID, c2, c4))
    fout.write("%s\n" % ( k2 ) )
    ID += 1
    hetePairs.add( (k1, k2, "nonsnp"+str(ID) ) )
ID = 0

for (k1, k2, c1, c2) in indelPair:
    fout.write(">kmer_snp%s_1_cov_%s\n" % (ID, c1))
    fout.write("%s\n" % ( k1 ) )
    fout.write(">kmer_snp%s_2_cov_%s\n" % (ID, c2))
    fout.write("%s\n" % ( k2 ) )
    ID += 1
    hetePairs.add( (k1, k2, "indel"+str(ID) ) )

fout.close()

left_index, right_index = kmercalling.build_left_right_kmer_index(uniqKmer)
extendPair = kmercalling.extend_pair(hetePairs, left_index, right_index, k)


f2 = "k_" + str(k)+"_extend_pair.snp"
with open(f2, "w") as fout:
    for (k1, k2, ID) in extendPair:
        fout.write(">kmer_%s_1\n" % (ID) )
        fout.write("%s\n" % ( k1 ) )
        fout.write(">kmer_%s_2\n" % (ID) )
        fout.write("%s\n" % ( k2 ) )


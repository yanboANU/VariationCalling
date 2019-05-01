#########################################################################
# File Name: variationCalling.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Mon 18 Mar 2019 13:59:49 AEDT
#########################################################################
#!/bin/bash
import argparse

from libprism.local import columns, contig
from libprism.local.prepare import clouds_from_refhap, merge_clouds, print_clouds
from libprism.local import tools
from math import log, exp
from pysam
import logging

#############################################################################

parser = argparse.ArgumentParser(description="Variation calling for Single-individual")

parser.add_argument('--bam', help='path to alignment bam file', required=True)
parser.add_argument('--ref', help='reference or scaffolds', required=True)
parser.add_argument('--a1', help='delete rate', required=True)
parser.add_argument('--a2', help='insert rate', required=True)
parser.add_argument('--obLen', help='clouds to chromosomes', required=True)
parser.add_argument('--cov', help='average coverage', required=True)

args = parser.parse_args()
logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

#############################################################################


with open(args.ref) as f:
    contigs = contig.read_Contig(f)

bamfile = pysam.AligmentFile(args.bam, "rb")

logging.info( "finish reading contigs and bam" )

a1 = float(args.a1)
a2 = float(args.a2) 
obLen = int(args.obLen) #defalt = 3
averageCov = int(args.cov)

DealLength = 5000000
Overlap = 150000
Divide = False
        
t = threading.Thread(target=loop, name='LoopTHread', args=(bamfile,contigs,a1,a2,obLen, averageCov))
t.start()
t.join()


def loop(bamfile, contigs, a1, a2, obLen, averageCov):

    for (contigName, contig) in contigs.items():
        print ("dealing ", contigName)
        if (contig._len < DealLength) or (Divide == False):
            deal_one_part(bamfile, contig, 0, contig._len, a1, a2, obLen, averageCov)
        else:
            t2 = threading.Thread(target=loop2, name='Loop2THread', args=(bamfile,contig,a1,a2,obLen, averageCov))
            t2.start()
            t2.join()
            
def loop2(bamfile, contig, a1, a2, obLen, averageCov):
    start = 0
    
    while start < contig._len:
        end = min(start+DealLength, contig._len)
        deal_one_part(bamfile, contig, start, end, a1, a2, obLen, averageCov)
        start = start + DealLength - Overlap

        
def deal_one_part(bamfile, contig, start, end, a1, a2, obLen, averageCov):
    time1 = time.clock()
    columns = column.init_Columns(bamfile, contig, False, start, end, a1, a2)
    time2 = time.clock()
    print ( "init columns one part running %s Seconds" % (time2 - time1) )
    p = phasing.Phasing(columns, contig, obLen)
    p._pre_Process(start, end, contig._len, averageCov)
    
    time3 = time.clock()
    print ( "phasing pre_process one part running %s Seconds" % (time3 - time2) )
    print ( "pre process finish" )
    
    p._phasing()
    time4 = time.clock()
    print ( "phasing one part running %s Seconds" % (time4 - time3) )

    print ( "finish phasing" )
    p._write_result2(start, end) 

    time5 = time.clock()
    print ( "writing one part running %s Seconds" % (time5 - time4) )    

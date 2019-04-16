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

#############################################################################

parser = argparse.ArgumentParser(description="Variation calling for Single-individual")

parser.add_argument('--bam', help='path to alignment bam file', required=True)
parser.add_argument('--ref', help='reference or scaffolds', required=True)
parser.add_argument('--a1', help='delete rate', required=True)
parser.add_argument('--a2', help='insert rate', required=True)
parser.add_argument('--obLen', help='clouds to chromosomes', required=True)
parser.add_argument('--cov', help='average coverage', required=True)

args = parser.parse_args()

#############################################################################

with open(args.ref) as f:
    contigs = contig.read_Contig(f)



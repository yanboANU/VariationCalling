#########################################################################
# File Name: test_run_tools.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Tue 16 Apr 2019 14:07:51 AEST
#########################################################################
#!/bin/bash

import os
import sys
import run_tools

bam_file = "/home/yulin/bio/HA/HapCUT2/pacbio/hg19/flag0flag16/5x/chr22.5x.bam"
ref_file = "/home/yulin/liyanbo/Data/reference/GRCh37_hg19/chr22.fa"
path_to_whatshap = "/home/yulin/bio/VariationCalling/simple-snp-caller.py"

run_tools.run_whatsHap(path_to_whatshap, ref_file, bam_file, "whatshap.vcf")

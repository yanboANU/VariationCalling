#########################################################################
# File Name: test_evaluate.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Tue 16 Apr 2019 14:41:33 AEST
#########################################################################
#!/bin/bash

import os
import sys
import error_rates


real_snp_file="/home/yanbo/bio/data/SNP_calling_result/NA12878/1000genome/22_snps"
high_snp_file="/home/yanbo/bio/data/SNP_calling_result/NA12878/High-confidence/22_snps"

error_rates.calc_TP(real_snp_file, high_snp_file, "whatshap.vcf")
error_rates.calc_TP(real_snp_file, high_snp_file, "chr22_snp_mutation_0_51304566")



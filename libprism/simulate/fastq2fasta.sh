#########################################################################
# File Name: fastq2fasta.sh
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Wed 03 Jul 2019 13:10:15 AEST
#########################################################################
#!/bin/bash

cat single_dat.fq | awk 'BEGIN{i=0}{ if(i%4==0 || i%4==1){ print $0} i+=1;}'

#########################################################################
# File Name: runkmercalling.sh
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Thu 04 Jul 2019 16:52:18 AEST
#########################################################################
#!/bin/bash
set -e
#cd $$3x/chr$1/
#mkdir kmercalling
#cd kmercalling

if [ $# != 5 ]; then
	echo "\$1:chrID \$2:k \$3:k-1 \$4:hom cov \$5:*fq"
	exit 
fi	

date
/home/yulin/software/dsk/build/bin/dsk  -file $5 -histo 1 -out chr$1_k$2 -kmer-size $2
/home/yulin/software/dsk/build/bin/dsk2ascii -file chr$1_k$2 -out chr$1_k$2.txt

/home/yulin/software/dsk/build/bin/dsk  -file $5 -histo 1 -out chr$1_k$3 -kmer-size $3
/home/yulin/software/dsk/build/bin/dsk2ascii -file chr$1_k$3 -out chr$1_k$3.txt

#findGSE:
/home/yulin/software/R-3.6.0/bin/Rscript /home/yulin/bio/VariationCalling/libprism/local/runfindGSE.r chr$1_k$2.histo $2 ./ $4 >hete.peak.k$2

left=`cat hete.peak.k$2 | grep "het_xfit_left" | awk '{print $6}'`
right=`cat hete.peak.k$2 | grep "het_xfit_right" | awk '{print $6}'`
echo $left
echo $right

python /home/yulin/bio/VariationCalling/variationCalling.py --t1 chr$1_k$2.txt --t2 chr$1_k$3.txt --c1 $left --c2 $right --k $2 >vc_k$2.log 

date


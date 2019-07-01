#########################################################################
# File Name: compareRealRair.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Thu 09 May 2019 11:51:38 AEST
#########################################################################
#!/bin/bash
import sys
import tools

def read_kmer(filename):
    kmer = []
    with open(filename, "r") as f:
        for line in f:
            words= line.strip().split()
            kmer.append(words)
    return kmer      



TNFile="uniq_sortedTN"
allFile="../30x_downsample/k" + sys.argv[1]  + "/chr22.sorted.adjust.txt"

list2 = read_kmer(TNFile)
list4 = read_kmer(allFile)

intersection=[]
index2, index4=0,0
len2, len4 = len(list2), len(list4)



while index2 < len2 and index4 < len4:
    if (list2[index2][0] < list4[index4][0]):
        print list2[index2][0], "not in all kmer( dsk )"
        index2 +=1
    elif list2[index2][0] > list4[index4][0]:   
        index4 +=1
    elif list2[index2][0] == list4[index4][0]:
        tools.print_list( list2[index2] )
        tools.print_list( list4[index4] )
        intersection.append( list4[index4] )
        index2 +=1
        index4 +=1
    else: 
        print "unknow"
        print list2[index2]
        print list4[index4]
        index2 +=1
        index4 +=1
        
print "TN kmer number, all kmer number, intersection pair kmer number"
print len2, len4, len(intersection)        
        


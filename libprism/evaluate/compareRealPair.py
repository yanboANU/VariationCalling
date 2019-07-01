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
    pair_kmer = []
    with open(filename, "r") as f:
        for line in f:
            words= line.strip().split()
            ''' 
            if len(words) == 2:
                pair_kmer.append((words[0], words[1]))
            if len(words) == 4:
                pair_kmer.append((words[0], words[1]))
            '''
            pair_kmer.append(words)
    return pair_kmer      



#realFile="chr22.real." + sys.argv[1] + "mer"
#pairFile="../30x_downsample/k" + sys.argv[1]  + "/chr22.pair.kmer"

realFile=sys.argv[1]
pairFile=sys.argv[2]


list2 = read_kmer(realFile)
list4 = read_kmer(pairFile)

intersection=[]
index2, index4=0,0
len2, len4 = len(list2), len(list4)
foutFP=open("FPPair", "w") #not real SNP, but found
foutTN=open("TNPair", "w") #real SNP, not found



while index2 < len2 and index4 < len4:
    if (list2[index2][0:2] < list4[index4][0:2]):
        foutTN.write( "%s %s\n" % (list2[index2][0], list2[index2][1] ) )
        index2 +=1
    elif list2[index2][0:2] > list4[index4][0:2]:   
        foutFP.write( "%s %s %s %s\n" % ( list4[index4][0], list4[index4][1], list4[index4][2], list4[index4][3] ) )
        index4 +=1
    elif list2[index2][0:2] == list4[index4][0:2]:
        #tools.print_list( list2[index2] )
        #tools.print_list( list4[index4] )
        intersection.append( list4[index4] )
        index2 +=1
        index4 +=1
    else: 
        print "unknow"
        print list2[index2]
        print list4[index4]
        index2 +=1
        index4 +=1
        
print "real SNP pair kmer number, find pair kmer number, intersection pair kmer number"
print len2, len4, len(intersection)        
        


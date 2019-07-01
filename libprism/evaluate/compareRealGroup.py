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
            pair_kmer.append(words)
    return pair_kmer      


def read_group_kmer(filename):
    pair_kmer = []
    pair_group_length = []
    state = 0
    with open(filename, "r") as f:
        for line in f:
            words= line.strip().split()
            if state == 0 and line.startswith("group"):
                pair_kmer.append(words[1:])
                state = 1
            elif state == 1:
                a = len(words)/2
                state = 2
            elif state == 2:
                b = len(words)/2
                state = 0
                pair_group_length.append( (a, b ) )

    return pair_kmer, pair_group_length      




#realFile="chr22.real.21mer"
#pairFile="chr22.group.kmer"

print "input: chr*.real.21mer chr*.group.kmer"
realFile=sys.argv[1]
pairFile=sys.argv[2]

real_pair_kmer = read_kmer(realFile)
pair_kmer, pair_group_length = read_group_kmer(pairFile)

intersection=[]
index2, index4=0, 0
len2, len4 = len(real_pair_kmer), len(pair_kmer)


foutFP = open("FPPair", "w") #not real SNP, but found
foutTP = open("TPPair", "w") #real SNP, not found
foutTN = open("TNPair", "w") #real SNP, not found

fout = open("GroupID2RefPos2", "w")
while index2 < len2 and index4 < len4:
    if (real_pair_kmer[index2][0:2] < pair_kmer[index4][1:3]):
        foutTN.write( "%s %s\n" % (real_pair_kmer[index2][0], real_pair_kmer[index2][1] ) )
        index2 +=1
    elif real_pair_kmer[index2][0:2] > pair_kmer[index4][1:3]:   
        foutFP.write( "%s %s %s %s " % ( pair_kmer[index4][1], pair_kmer[index4][2], pair_kmer[index4][3], pair_kmer[index4][4] ) )
        foutFP.write( "%s %s\n" % ( pair_group_length[index4][0], pair_group_length[index4][1] ) )
        index4 +=1
    elif real_pair_kmer[index2][0:2] == pair_kmer[index4][1:3]:
        #tools.print_list( real_pair_kmer[index2] )
        #tools.print_list( pair_kmer[index4] )
        intersection.append( pair_kmer[index4] )
        foutTP.write( "%s %s %s %s %s " % ( pair_kmer[index4][1], pair_kmer[index4][2], real_pair_kmer[index2][2], pair_kmer[index4][3], pair_kmer[index4][4] ) )
        foutTP.write( "%s %s\n" % ( pair_group_length[index4][0], pair_group_length[index4][1] ) )
        fout.write("%s %s\n" % (pair_kmer[index4][0], real_pair_kmer[index2][2] ) )
        index2 +=1
        index4 +=1
    else: 
        print "unknow"
        print real_pair_kmer[index2]
        print pair_kmer[index4]
        index2 +=1
        index4 +=1
fout.close()
foutTP.close()
foutFP.close()   
foutTN.close()
TPnum = len(intersection)
print "real SNP pair kmer number, find group kmer number, intersection pair kmer number, accuracy, sensitive"
print len2, len4, TPnum, round(float(TPnum)/len4,3), round(float(TPnum)/len2,3)
        


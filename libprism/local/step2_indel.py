#########################################################################
# File Name: step2.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Wed 08 May 2019 11:05:49 AEST
#########################################################################
#!/bin/bash
# find pair k-mer(21-mer) 

import tools
import sys

def find_pair_indel_kmer(kmerFile, k_1merFile, output_filename):    
    kmers = {}
    count = 0
    with open(kmerFile, "r") as f:
        for line in f:
            words = line.strip().split()
            l = len(words[0])
            newtemp = tools.reverse(words[0]) # these two lines 
            assert newtemp >= words[0]        # can be deleted  

            mid = int(l/2)
            leftHalf = words[0][:mid]
            rightHalf = words[0][mid+1:]

            if (tools.hamming_distance(rightHalf, words[0][mid : -1 ] ) <= 1 or # mutation => delete 
                    tools.hamming_distance(leftHalf, words[0][1:mid+1]) <= 1 ):
                continue

            #if words[0][mid] == words[0][mid-1] and words[0][mid] == words[0][mid+1] :  # only keep those pair, delete at first position or last position
            #if words[0][mid] == words[0][mid-1]:  # only keep those pair, delete at first position or last position
                #continue

            key = leftHalf + rightHalf

            Rkey = tools.reverse(key)
            assert key <= Rkey
            if key not in kmers:
                kmers[key] = []
            kmers[key].append( (words[0], words[1]) )
            count += 1

    print ("unique kmer number", count)
    print ("unique k_1mer in kmers", len(kmers) )

    k_1mers = {}
    count = 0
    with open(k_1merFile, "r") as f:
        for line in f:
            words = line.strip().split()
            l = len(words[0])
            key = words[0]
            newtemp = tools.reverse(key) # these two lines 
            assert newtemp >= key       # can be deleted  
            assert key not in k_1mers 
            k_1mers[key] =  words[1] # words[1] is coverage
            count += 1
    print ("unique k_1mer number", count)

    pairKmers = []
    for key in k_1mers:
        if key in kmers:
            if len(kmers[key]) == 1:
                #print key, k_1mers[key], kmers[key]
                kmer, c = kmers[key][0]
                if kmer[:-1] == key or kmer[1:] == key: # remove head or tail, they are same
                    continue
                pairKmers.append( (kmer, key, c, k_1mers[key]) )
            else:
                print ("more than 1", key, k_1mers[key], kmers[key])

    
    fout = open(output_filename, "w")
    sortedKmers = sorted(pairKmers)
    countU, countD, count = 0, 0, 0
    for (k1, k2, c1, c2) in sortedKmers:
        #print (mid)
        if k1[mid] == k1[mid+1] or k1[mid] == k1[mid-1]:
            countD += 1
        else:
            countU += 1
        count += 1 
        fout.write("%s %s %s %s\n" % ( k1, k2, c1, c2 ) )
    fout.close()       
    print (countU, countD, count)
    


kmerFile = sys.argv[1] # k.unique.kmer
k_1merFile = sys.argv[2] # k_1.unique.kmer
output_filename= sys.argv[3] # ".pair.kmer"
find_pair_indel_kmer(kmerFile, k_1merFile, output_filename)

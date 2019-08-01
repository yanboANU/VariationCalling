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

def find_pair_kmer(input_filename, output_filename):    
    m={}
    with open(input_filename, "r") as f:
        count, count2 = 0,0

        for line in f:
            #if count >100:
                #break
            words = line.strip().split()
            l = len(words[0])
            newtemp = tools.reverse(words[0]) # these two lines 
            assert newtemp >= words[0]        # can be deleted  
            key= words[0][:int(l/2)] + words[0][int(l/2)+1:]
            if key not in m:
                m[key] = []
            m[key].append( (words[0], words[1]) )
            count += 1

    print ("unique kmer number", count)        
    print ("not small kmer (compare reverse kmer)", count2)
    print ("total number possible pair kmer", len(m))    

    pairKmers = []
    count1, countLarge2 =0, 0
    for key in m:
        if len(m[key]) == 1:
            count1 += 1
        elif len(m[key]) == 2:
            k1 = m[key][0][0]
            k2 = m[key][1][0]
            newk1 = tools.reverse(k1)             # these three line 
            newk2 = tools.reverse(k2)             # can be deleted
            assert k1 <= newk1 and k2 <= newk2
            #min_k = min(k1, k2)
            #min_new = min(newk1, newk2)
            #if min_new < min_k:
                #k1 = newk1
                #k2 = newk2
            #sum_coverage = int(m[key][0][1]) + int (m[key][1][1])
            if k1 < k2:
                pairKmers.append( (k1, k2, m[key][0][1], m[key][1][1]) )
            else:
                pairKmers.append( (k2, k1, m[key][1][1], m[key][0][1]) )
        else:
            countLarge2 += 1

    print ("kmer cannot find pair number", count1)       
    print ("more than one mutation in middle", countLarge2)

    fout = open(output_filename, "w")
    sortedKmers = sorted(pairKmers)
    for (k1, k2, c1, c2) in sortedKmers:      
        fout.write("%s %s %s %s\n" % ( k1, k2, c1, c2 ) )
    fout.close()        



input_filename= sys.argv[1] + ".unique.kmer"
output_filename= sys.argv[1] + ".pair.kmer"
find_pair_kmer(input_filename, output_filename)

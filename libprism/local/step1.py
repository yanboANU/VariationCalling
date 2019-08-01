#########################################################################
# File Name: step0.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Wed 08 May 2019 11:05:49 AEST
#########################################################################
#!/bin/bash
import sys
import tools

# dsk not always use smaller kmer (compare with reverse kmer)
def pick_smaller_unique_kmer(input_filename, low, high, output_filename):
    fout = open(output_filename, "w")
    with open(input_filename, "r") as f:
        for line in f:
            words = line.strip().split()
            coverage = int(words[1])
            if coverage < low or coverage > high:
                continue 
            kmer = words[0]
            newkmer = tools.reverse(kmer)
            if kmer > newkmer:    
                kmer = newkmer
            fout.write("%s %s\n" % (kmer, words[1]))
        fout.close()        

input_filename = sys.argv[1] + ".txt"
output_filename = sys.argv[1] + ".unique.kmer"
low = int(sys.argv[2])
high = int(sys.argv[3])
print "input:", input_filename, "coverage range", low, high
pick_smaller_unique_kmer(input_filename, low, high, output_filename)

'''
def pick_smaller_kmer(input_filename, output_filename):
    fout = open(output_filename, "w")
    with open(input_filename, "r") as f:
        for line in f:
            words = line.strip().split()
            kmer = words[0]
            newkmer = tools.reverse(kmer)
            if kmer > newkmer:    
                kmer = newkmer
            fout.write("%s %s\n" % (kmer, words[1]))
        fout.close()        
input_filename = sys.argv[1] + ".txt"
output_filename = sys.argv[1] + ".adjust.txt"
pick_smaller_kmer(input_filename, output_filename)
'''

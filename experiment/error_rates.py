#########################################################################
# File Name: error_rate.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Tue 16 Apr 2019 14:40:27 AEST
#########################################################################
#!/bin/bash

import os
import sys


        


def read_vcf2(filename): 
    f = open(filename, "r")
    insert1s, delete1s, snps = {}, {} , {}
    for line in f:
        if line.startswith('#'):
            continue
        words = line.split()
        chrID = words[0]
        if chrID not in snps:
            insert1s[chrID] , delete1s[chrID] , snps[chrID] = [], [], []
        s1 = words[3].strip()
        s2 = words[4].strip()

        if s1 == '.' or s2 == '.':
            continue
        homo = words[9].split(':')[0]
        s1Len = len(s1)
        s2Len = len(s2)
        
        if s1Len == 1 and s2Len == 1:
            snps[chrID].append((words[1], s1, s2, homo))
            continue
                
        if s1Len == 1 and s2Len >= 3:
            sub = s2.split(',')
            subLen = len(sub)
            if subLen >= 2:
                i=0
                flag = True
                for i in range(subLen):
                    if len(sub[i]) != 1:
                        flag = False
                if flag == True:
                    snps[chrID].append((words[1], s1, s2, homo))
                    continue

        if s1Len == 1 and s2Len >= 2 and s2[0] == s1:
            if s2Len == 3 and s2.count(',')>=1:
                print s2
                continue
            insert1s[chrID].append((words[1], s1, s2, homo))
            continue        
           
        if s1Len >= 2 and s2Len == 1 and s1[0] == s2:
            delete1s[chrID].append((words[1], s1, s2, homo))
            continue
    f.close()       
    return insert1s, delete1s, snps    
 

def read_vcf(filename): # only include snp
    pos = []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith('#'):
                continue
            words = line.split()
            if words[1].isdigit():
                pos.append(int(words[1]))
    return pos


def read_int_list(filename): # read_mutation_list
    f = open(filename, "r")
    pos = []  
    for line in f:  
        words = line.split()
        #if len(words) == 1 or len(words) == 3 or len(words) == 5:
        if words[0].isdigit():
            pos.append(int(words[0]))
    return pos

#real_snp_file=/home/yanbo/bio/data/SNP_calling_result/NA12878/1000genome/22_snps
#high_snp_file=/home/yanbo/bio/data/SNP_calling_result/NA12878/High-confidence/22_snps

def calc_TP(real_snp_file, high_snp_file, find_snp_file):
    
    
    # position in ref
    real_snp_pos = read_int_list(real_snp_file) 
    # position in high
    high_snp_pos = read_int_list(high_snp_file)
    # position in find
    find_snp_pos = read_vcf(find_snp_file)
 
    print ("number of real snp:", len(real_snp_pos))
    print ("number of high snp:", len(high_snp_pos))
    print ("number of find snp:", len(find_snp_pos))
     
    pr = set(find_snp_pos).intersection(set(real_snp_pos))
    ph = set(find_snp_pos).intersection(set(high_snp_pos))
    hr = set(real_snp_pos).intersection(set(high_snp_pos))
    phr = pr.intersection(set(high_snp_pos))
    print ("pr number:", len(pr))
    print ("ph number:", len(ph))
    print ("hr number:", len(hr))
    print ("phr number:", len(phr))
    #print ("FP:", sorted(set(find_snp_pos) - set(real_snp_pos)))
    #print ("TN number:", len(set(real_snp_pos)- set(find_snp_pos)))
    #print ("TN:", sorted(set(real_snp_pos)- set(find_snp_pos)))




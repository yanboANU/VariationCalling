#########################################################################
# File Name: plot_gap_distribution.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Fri 22 Mar 2019 15:13:10 AEDT
#########################################################################
#!/bin/bash
import matplotlib
#import numpy as np
import matplotlib.pyplot as plt
import sys
import math

#use one vector, real histogram

def read_data(filename, label):
    cov, covMinus = [], []
    covSum, covDiv = [], []
    gLen, gMinus = [], []
    gSum, gDiv = [], []
    with open(filename) as f:
        for line in f:
            words = line.strip().split()
            if label == "FP":
                a,b = float(words[2]),float(words[3])
            if label == "TP":
                a,b = float(words[2]),float(words[3])
            cov.append(a)
            cov.append(b)
            covSum.append(a+b)
            covMinus.append(a-b)
            covDiv.append( round(a/b, 4) )
            '''
            gLen.append(c)
            gLen.append(d)
            if c+d >102:
                print (line, c,d, c+d)
            assert c+d<=102
            gSum.append(c+d)
            gMinus.append(c-d)
            gDiv.append( round(c/d, 4) )
            '''
    return cov, covSum, covMinus, covDiv #, gLen, gSum, gMinus, gDiv



TPfilename=sys.argv[1] # TPPair

TPcov, TPcovSum, TPcovMinus, TPcovDiv = read_data(TPfilename, "TP")
#TPcov, TPcovSum, TPcovMinus, TPcovDiv, TPgLen, TPgSum, TPgMinus, TPgDiv = read_data(TPfilename, "TP")
FPfilename=sys.argv[2] # FPPair
#FPcov, FPcovSum, FPcovMinus, FPcovDiv, FPgLen, FPgSum, FPgMinus, FPgDiv = read_data(FPfilename, "FP")
FPcov, FPcovSum, FPcovMinus, FPcovDiv = read_data(FPfilename, "FP")

plt.hist(TPcovMinus, 50, histtype='step', label='TP')
plt.hist(FPcovMinus, 50, histtype='step', label='FP')
plt.legend()
plt.xlabel('Pair Kmer Coverage Minus', fontsize=18)
plt.ylabel('Count', fontsize=18)
plt.savefig("cov_minus_dis.png")

plt.clf()

plt.hist(TPcovDiv, 10, histtype='step', label='TP')
plt.hist(FPcovDiv, 10, histtype='step', label='FP')
plt.legend()
plt.xlabel('Pair Kmer Coverage Div', fontsize=18)
plt.ylabel('Count', fontsize=18)
plt.savefig("cov_div_dis.png")

'''
plt.clf()
plt.hist(TPgMinus, 10, histtype='step', label='TP')
plt.hist(FPgMinus, 10, histtype='step', label='FP')
plt.legend()
plt.xlabel('Pair kmer Group size Minus', fontsize=18)
plt.ylabel('Count', fontsize=18)
plt.savefig("group_size_minus_dis.png")

plt.clf()
#print (TPgDiv)
#print (FPgDiv)
plt.hist(TPgDiv, 40, histtype='step', label='TP')
plt.hist(FPgDiv, 40, histtype='step', label='FP')
plt.legend()
plt.xlim(0.002,5)
plt.xlabel('Pair Kmer group size Div', fontsize=18)
plt.ylabel('Count', fontsize=18)
plt.savefig("group_size_div_dis.png")
'''
plt.clf()
plt.hist(TPcov, 50, histtype='step', label='TP')
plt.hist(FPcov, 50, histtype='step', label='FP')
plt.legend()
plt.xlabel('21-mer Coverage', fontsize=18)
plt.ylabel('Count', fontsize=18)
plt.savefig("cov_dis.png")

plt.clf()
plt.hist(TPcovSum, 50, histtype='step', label='TP')
plt.hist(FPcovSum, 50, histtype='step', label='FP')
plt.legend()
plt.xlabel('Pair Kmer Coverage Sum', fontsize=18)
plt.ylabel('Count', fontsize=18)
plt.savefig("cov_sum_dis.png")

'''
plt.clf()
plt.hist(TPgLen, 21, histtype='step', label='TP', log=True)
plt.hist(FPgLen, 21, histtype='step', label='FP', log=True)
plt.ylim(1, 100000)
plt.legend()
plt.xlabel('Group size', fontsize=18)
plt.ylabel('Count', fontsize=18)
plt.savefig("group_size_dis.png")

plt.clf()
plt.hist(TPgSum, 20, histtype='step', label='TP', log=True)
plt.hist(FPgSum, 20, histtype='step', label='FP', log=True)
plt.legend()
plt.xlabel('Pair Kmer group size Sum', fontsize=18)
plt.ylabel('Count', fontsize=18)
plt.savefig("group_size_sum_dis.png")
'''

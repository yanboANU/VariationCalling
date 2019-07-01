#########################################################################
# File Name: tools.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Thu 09 May 2019 13:12:45 AEST
#########################################################################
#!/bin/bash

import sys
import os
#import subprocess

m={}
m['A'] = 'T'
m['T'] = 'A'
m['C'] = 'G'
m['G'] = 'C'

def reverse(s): # AATCG
    news = ""
    for ele in s:
       news = m[ele] + news #CGATT
    return news   

def reverse_ward(ward):
    if ward == 'f':
        return 'b'
    elif ward == 'b':
        return 'f'
    else:
        print "reverse forward or backward error"
        sys.exit()

'''
s="ATCG"
print s
print reverse(s)
'''
def print_list(l):
    for ele in l:
        print ele,
    print ""    

def file_lines(filename):
    command = "wc -l " + filename + " >file_lines"
    os.system(command)
    with open("file_lines", "r") as f:
        for line in f:
            words = line.split()
            return int(words[0])


def hamming_distance(s1, s2):
    count = 0
    lenS = len(s1)
    #print s1, s2
    assert lenS == len(s2)
    for i in range(lenS):
        if s1[i] == -1 or s2[i] == -1:
            continue
        if s1[i] != s2[i]:
            count +=1
    return count

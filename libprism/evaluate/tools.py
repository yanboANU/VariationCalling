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

def reverse(s): # ATCG
    news = ""
    for ele in s:
       news = m[ele] + news
    return news   

'''
s="ATCG"
print s
print reverse(s)
'''
def print_list(l):
    for ele in l:
        print (ele, )
    print ("")

def file_lines(filename):
    command = "wc -l " + filename + " >file_lines"
    os.system(command)
    with open("file_lines", "r") as f:
        for line in f:
            words = line.split()
            return int(words[0])

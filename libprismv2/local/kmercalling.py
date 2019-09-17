#########################################################################
# File Name: kmercalling.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Fri 09 Aug 2019 11:40:22 AEST
#########################################################################
#!/bin/bash

import os
import sys
import time
from libprismv2.local import tools
#from tools import *

def pick_smaller_unique_kmer(input_filename, low, high):
    uniqKmer = []
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
            uniqKmer.append( (kmer, coverage) )
    return uniqKmer

# hamming distance = 1
def find_snp_pair_kmer(uniqKmers, k, left_index, right_index):    
    m={}
    mid = int(k/2)
    print ("before build mapK")
    for (kmer, cov) in uniqKmers:
        #Rkmer = tools.reverse(kmer) # these two lines 
        #assert Rkmer >= kmer        # can be deleted  
        key= kmer[:mid] + kmer[mid+1:]
        if key not in m:
            m[key] = []
        m[key].append( (kmer, cov) )
    print ("after build mapK")
    print ("unique kmer number", len(uniqKmers))        
    print ("total number possible pair kmer", len(m))    
    pairKmers = []
    count1, countLarge2 =0, 0
    for key in m:
        mKeyLen = len(m[key])
        if mKeyLen == 1:
            count1 += 1
        elif mKeyLen == 2:
            k1 = m[key][0][0]
            k2 = m[key][1][0]
            if k1[1:] == k2[:-1] or k1[:-1]==k2[1:]:
                print (k1, k2)
                continue
            # if you want to find more, don't use this
            # use this, find less but high accuracy
            ek1, ek2, flag = extend_one_pair(k1, k2, left_index, right_index,k, 0)
            if flag == False:
                continue
            if k1 < k2:
                pairKmers.append( (k1, k2, ek1, ek2, m[key][0][1], m[key][1][1]) )
            else:
                pairKmers.append( (k2, k1, ek2, ek1, m[key][1][1], m[key][0][1]) )
        else:
            countLarge2 += 1
    
    print ("kmer cannot find pair number", count1)       
    print ("more than one mutation in middle", countLarge2) 
    fout = open("snp_pair", "w")
    sortedKmers = sorted(pairKmers)
    ID = 1
    for (k1, k2, ek1, ek2, c1, c2) in sortedKmers:   
        fout.write("%s %s %s %s %s %s\n" % (k1, k2, ek1, ek2, c1, c2))
    fout.close()        
    return m, pairKmers
    

'''
print ("before filter mapK")
uniqMapK = {}
for key in mapK:
    for (kmer, cov) in mapK[key]:
        leftHalf = kmer[:mid]
        rightHalf = kmer[mid+1:]
        if (tools.hamming_distance(rightHalf, kmer[mid : -1 ] ) <= 1 or # mutation => delete 
                tools.hamming_distance(leftHalf, kmer[1:mid+1]) <= 1 ):
            mapK[key].remove((kmer, cov))
    if len(mapK[key]) == 1:
        uniqMapK[key] = mapK[key]

print ("after filter mapK")
'''
def find_indel_pair_kmer(uniqKmers, uniqK_1mers, k, left_index, right_index, threshold): #, highRepeat):    
    mid = int(k/2)
    mapK = {}
    for (temp, cov) in uniqKmers:
        leftHalf = temp[:mid]
        rightHalf = temp[mid+1:]
        
        if (tools.hamming_distance(rightHalf, temp[mid : -1 ] ) <= 2 or # mutation => delete 
                tools.hamming_distance(leftHalf, temp[1:mid+1]) <= 2 ):
            continue
        #if temp in highRepeat:
        #    continue
        #GC = (temp.count('C')+temp.count('G'))/float(k)
        #if GC < 0.05 or GC > 0.85:
        #    continue
        key = leftHalf + rightHalf
        #Rkey = tools.reverse(key)
        #assert key <= Rkey
        if key not in mapK:
            mapK[key] = []
        mapK[key].append( (temp, cov) )
    mapK_1 = {} 
    print ( "uniq K-1 mer size", len(uniqK_1mers) )
    for (key, cov) in uniqK_1mers:        
        mapK_1[key] = cov
    indelPair = []
    print ( "map K-1 size", len(mapK_1) )
    print ( "check", k/4)
    for key in mapK_1:
        if key in mapK and len( mapK[key] ) == 1:
            kmer, c = mapK[key][0]
            ek1, ek2, flag = extend_one_pair(kmer, key, left_index, right_index,k, threshold)
            if flag == False:
                continue
            indelPair.append( (kmer, key, ek1, ek2, c, mapK_1[key]) )
     
    fout = open("indel_pair", "w")
    sortedKmers = sorted(indelPair)
    for (k1, k2, ek1, ek2, c1, c2) in sortedKmers:
        fout.write("%s %s %s %s %s %s\n" % (k1, k2, ek1, ek2, c1, c2))
    fout.close()       
    return indelPair

def update_map_merge(left_i, left_j, key1, key2, mapMerge):
    
    k1, cov1 = left_i
    k2, cov2 = left_j
    '''
    if key2 < key1:
        key1, key2 = key2, key1
        k1, k2 = k2, k1
        cov1, cov2 = cov2, cov1
    if (key1, key2) not in mapMerge:
        mapMerge[ (key1, key2) ] = []  
    mapMerge[ (key1, key2) ].append( (k1, k2, cov1, cov2) )
    '''
    minkey1, minkey2 = tools.get_smaller_pair_kmer(key1, key2)
    if (minkey1, minkey2) not in mapMerge:
        mapMerge[ (minkey1, minkey2) ] = []   
    Rk1, Rk2 = tools.reverse(k1), tools.reverse(k2)
    if k1.count(minkey1) == 1 and k2.count(minkey2) == 1:
        mapMerge[ (minkey1, minkey2) ].append( (k1, k2, cov1, cov2) )
    elif k1.count(minkey2) == 1 and k2.count(minkey1) == 1: 
        mapMerge[ (minkey1, minkey2) ].append( (k2, k1, cov2, cov1) )
    elif Rk1.count(minkey1) == 1 and Rk2.count(minkey2) == 1:   
        mapMerge[ (minkey1, minkey2) ].append( (Rk1, Rk2, cov1, cov2) )
    elif Rk1.count(minkey2) == 1 and Rk2.count(minkey1) == 1:
        mapMerge[ (minkey1, minkey2) ].append( (Rk2, Rk1, cov2, cov1) )
    else:
        print ("something wrong 1")
        sys.exit()
    
    return     

def build_map_merge(left, k):
    mapMerge= {}
    mid = int(k/2)
    hisMap = {}
    for key in left:
        groupSize = len(left[key])
        for i in range(0, groupSize-1):
            for j in range(i+1, groupSize):
                k1, cov1 = left[key][i]
                k2, cov2 = left[key][j]
                if k1[mid] == k2[mid]: #or (k1 in mappedKmer) or (k2 in mappedKmer):
                    continue
                dis = 0
                cnt = mid + 1
                diffPos = 0
                while cnt < k:
                    if k1[cnt] != k2[cnt]:
                        dis += 1
                        diffPos = cnt
                    cnt += 1    
                    if dis >= 2:
                        break
                if cnt == k and dis == 1:
                    key1, key2 = k1[ diffPos - mid : ] , k2[ diffPos - mid : ]
                    mink1, mink2 = tools.get_smaller_pair_kmer(k1, k2)
                    #fout.write("%s %s %s %s\n" % (mink1, mink2, cov1, cov2) )
                    #candidateNonPair.append((mink1, mink2))
                    if mink1 not in hisMap:
                        hisMap[mink1] = 0
                    if mink2 not in hisMap:
                        hisMap[mink2] = 0
                    hisMap[mink1] += 1
                    hisMap[mink2] += 1
                    update_map_merge(left[key][i], left[key][j], key1, key2, mapMerge)
                    #mappedKmer.add(k1)
                    #mappedKmer.add(k2)
                    #break #3 lines add 22 Aug. a kmer only allow one kmer hamming distance equal to 2
    
    # for debug, can be deleted            
    highRepeat = set()
    for key in hisMap:
        if hisMap[key] > 2:
            highRepeat.add(key)
    print ("high Repeat kmer number", len(highRepeat) )
    #fout.close()                
    return mapMerge, highRepeat

def merge_pair(mapMerge, highRepeat, k, left_index, right_index):
    pairSet, nonPair = set(), []
    for (key1, key2) in mapMerge:
        mlen = len( mapMerge [ (key1, key2) ] )
        if mlen > 2 or mlen==1: # one key only allow a pair
            continue
        Left1, Left2, covL1, covL2 = mapMerge[ (key1, key2) ][0]
        Right1, Right2, covR1, covR2 = mapMerge[ (key1, key2) ][1]
        l = len(key1)

        Lmin1 = min(Left1, tools.reverse(Left1))
        Lmin2 = min(Left2, tools.reverse(Left2))
        Rmin1 = min(Right1, tools.reverse(Right1))
        Rmin2 = min(Right2, tools.reverse(Right2))
        if (Lmin1 in highRepeat or Lmin2 in highRepeat or
                Rmin1 in highRepeat or Rmin2 in highRepeat):
            continue
        if covL1 > 1.5*covR1 or 1.5*covL1 < covR1:
            continue
        if covL2 > 1.5*covR2 or 1.5*covL2 < covR2:
            continue
        if Right1[-l:] == Left1[:l] and Right2[-l:] == Left2[:l]:
            merge1 = Right1 + Left1[l:]
            merge2 = Right2 + Left2[l:]
        elif Left1[-l:] == Right1[:l] and Left2[-l:] == Right2[:l]:
            merge1 = Left1 + Right1[l:]
            merge2 = Left2 + Right2[l:]
        else:
            #print ("one side")
            continue
         
        #TTTTTTTTTTTTTTTCAAAAAAAAAAAAAAAA # also useful SNP and indel 
        #TTTTTTTTTTTTTTTTCAAAAAAAAAAAAAAA 
        if merge1[1:] == merge2[:-1] or merge1[:-1]==merge2[1:]:
            print (merge1, merge2)
            continue
        small1, small2 = tools.get_smaller_pair_kmer(merge1, merge2)
        ek1, ek2, flag = extend_one_pair(small1, small2, left_index, right_index, k, 0)
        if flag == False:
            continue
        if (small1, small2) not in pairSet:
            pairSet.add( (small1, small2) )
            nonPair.append( (small1, small2, ek1, ek2, covL1, covL2, covR1, covR2) )
    return nonPair        


def find_non_pair_kmer(uniqKmer, k, left_index, right_index):

    mid = int(k/2)
    left = {}
    for (kmer, cov) in uniqKmer:
        leftKey = kmer[:mid] 
        if leftKey not in left:
            left[ leftKey ] = []  
        left[leftKey].append( (kmer, cov) )
        Rkmer = tools.reverse(kmer)
        leftKey = Rkmer[:mid] 
        if leftKey not in left:
            left[ leftKey ] = []  
        left[leftKey].append( (Rkmer,cov) )
    
    #build map: overlap is key
    print ("left size", len(left) )
    mapMerge, highRepeat = build_map_merge(left, k)
    print ("map Merge size", len(mapMerge) )
    nonPair = merge_pair(mapMerge, highRepeat, k, left_index, right_index) 

    print ("non pair size", len(nonPair) )
    fout = open("non_pair", "w")
    sortedNon = sorted(nonPair)
    print ("non pair size", len(sortedNon) )
    for (k1, k2, ek1, ek2, c1, c2, c3, c4) in sortedNon:     
        fout.write("%s %s %s %s %s %s %s %s\n" % (k1,k2,ek1,ek2,c1,c2,c3,c4) )
    fout.close()        
    
    return nonPair, highRepeat


def build_left_right_kmer_index(uniqKmer):
    
    left_index , right_index = {}, {}
    count = 0
    for (kmer, cov) in uniqKmer:
        key = kmer[:-1]
        if key not in left_index:
            left_index[key] = list()
        left_index[key].append(kmer[-1])
        key = kmer[1:]
        if key not in right_index:
            right_index[key] = list()
        right_index[key].append(kmer[0])

    left_uniq_index , right_uniq_index = {}, {}
    for key in left_index:
        if len(left_index[key]) == 1:
            left_uniq_index[key] = left_index[key]

    for key in right_index:
        if len(right_index[key]) == 1:
            right_uniq_index[key] = right_index[key]
    return left_uniq_index, right_uniq_index


def extend_to_left(h1, left_index, right_index, k):
    
    mid = int(k/2)
    key = h1[ : (k-1)]
    Rkey = tools.reverse(key)
    temp, Rtemp = h1, tools.reverse(h1)
    add, Radd = "", ""
    for i in range(0, mid):
        flag, flagR = False, False
        if key in right_index :
            temp = right_index[key][0] + temp
            Rtemp = tools.reverse(temp)
            flag = True

        if Rkey in left_index:
            Rtemp = Rtemp + left_index[Rkey][0]
            temp = tools.reverse(Rtemp)
            flagR = True

        if flag == True and flagR == False:
            add = right_index[key][0] + add
            Radd = tools.reverse(add) 
            key = temp[: (k-1) ]
            Rkey = tools.reverse(key)
        elif flag == False and flagR == True: 
            Radd = Radd + left_index[Rkey][0]
            add = tools.reverse(Radd)
            Rkey = Rtemp[-(k-1):]
            key = tools.reverse(Rkey)
        elif flag == flagR:
            if flag == True:
                temp = temp[1:]
                Rtemp = Rtemp[:-1]
            break
    return temp, add

def extend_to_right(h1, left_index, right_index, k):

    mid = int(k/2)
    key = h1[-(k-1):]
    Rkey = tools.reverse(key)
    temp, Rtemp = h1, tools.reverse(h1)
    add, Radd = "", ""
    for i in range(0, mid):
        flag, flagR = False, False
        if key in left_index:
            temp = temp + left_index[key][0]
            Rtemp = tools.reverse(temp)
            flag = True

        if Rkey in right_index:
            Rtemp = right_index[Rkey][0] + Rtemp
            temp = tools.reverse(Rtemp)
            flagR = True

        if flag == True and flagR == False:
            add = add + left_index[key][0]
            Radd = tools.reverse(add)
            key = temp[-(k-1):]
            Rkey = tools.reverse(key)
        elif flag == False and flagR == True: 
            Radd = right_index[Rkey][0] + Radd
            add = tools.reverse(Radd)
            Rkey = Rtemp[:(k-1)]
            key = tools.reverse(Rkey)
        elif flag == flagR:
            if flag == True:
                temp = temp[:-1]
                Rtemp = Rtemp[1:]
            break
    return temp, add
   
def extend_one_pair(h1,h2,left_index, right_index,k, threshold):

    flag = True 
    temp1, add1 = extend_to_left(h1, left_index, right_index, k)
    temp2, add2 = extend_to_left(h2, left_index, right_index, k)
    minL = min (len(add1), len(add2) )
    # some TP add content shift one position equal 
    if minL!= 0 and add1[-minL:] != add2[-minL:]:
        flag = False    
    if minL == int(k/2) and tools.hamming_distance(add1, add2) == 1:
    #if minL == int(k/2) and tools.min_edit_distance(add1, add2) <= 2:
        flag = True    
    ekmer1, add1R = extend_to_right(temp1, left_index, right_index, k)
    ekmer2, add2R = extend_to_right(temp2, left_index, right_index, k) 
    minR = min (len(add1R), len(add2R) )
    if minR!=0 and add1R[0:minR] != add2R[0:minR]:
        flag = False
    # the distance of two snps larger than k/2 smaller than k    
    if minR == int(k/2) and tools.hamming_distance(add1R, add2R) == 1:
    #if minR == int(k/2) and tools.min_edit_distance(add1R, add2R) <= 2:
        flag = True
    if max(minL, minR) <= 2: 
        flag = False
    if min(minL, minR) <= threshold:
        flag = False
    return ekmer1, ekmer2, flag

def extend_pair(hetePairs, left_index, right_index, k):
    extendPair, extendSet = [], set()
    for (h1, h2, ID) in hetePairs:
        ekmer1, ekmer2 = extend_one_pair(h1, h2, left_index, right_index,k)
        small1, small2 = tools.get_smaller_pair_kmer(ekmer1, ekmer2)
        if (small1, small2) not in extendSet:
            extendSet.add( (small1, small2) )
            extendPair.append( (small1, small2, ID) )
    return extendPair   

def check_unique_next(kmer, left_index, right_index):
    Rkmer = tools.reverse(kmer)
    if Rkmer < kmer:
        kmer = Rkmer
    leftKey = kmer[1:]
    rightKey = kmer[:-1]
    #if leftKey in left_index and rightKey in right_index:

    if (leftKey in left_index) or (rightKey in right_index):
        #assert len(left_index[leftKey]) == 1
        #assert len(right_index[rightKey]) == 1
        return True

    return False

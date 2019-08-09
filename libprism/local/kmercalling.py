#########################################################################
# File Name: kmercalling.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Fri 09 Aug 2019 11:40:22 AEST
#########################################################################
#!/bin/bash

import os
import sys
import tools
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
def find_snp_pair_kmer(uniqKmers, k):    
    m={}
    mid = int(k/2)
    print ("before build mapK")
    for (kmer, cov) in uniqKmers:
        Rkmer = tools.reverse(kmer) # these two lines 
        assert Rkmer >= kmer        # can be deleted  
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
        if len(m[key]) == 1:
            count1 += 1
        elif len(m[key]) == 2:
            k1 = m[key][0][0]
            k2 = m[key][1][0]
            if k1 < k2:
                pairKmers.append( (k1, k2, m[key][0][1], m[key][1][1]) )
            else:
                pairKmers.append( (k2, k1, m[key][1][1], m[key][0][1]) )
        else:
            countLarge2 += 1

    print ("kmer cannot find pair number", count1)       
    print ("more than one mutation in middle", countLarge2)
    '''
    fout = open("snp_pair", "w")
    sortedKmers = sorted(pairKmers)
    ID = 1
    for (k1, k2, c1, c2) in sortedKmers:     
        fout.write(">kmer_snp%s_1_cov_%s\n" % (ID, c1))
        fout.write("%s\n" % ( k1 ) )
        fout.write(">kmer_snp%s_2_cov_%s\n" % (ID, c2))
        fout.write("%s\n" % ( k2 ) )
        ID += 1
    fout.close()        
    '''
    return m, pairKmers
    

def find_indel_pair_kmer(mapK, uniqK_1mers, k):    
    mid = int(k/2)
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
    mapK_1 = {}

    print ( "uniq K-1 mer size", len(uniqK_1mers) )
    for (key, cov) in uniqK_1mers:        
        #assert key not in uniqK_1mers 
        mapK_1[key] = cov

    indelPair = []

    print ( "map K-1 size", len(mapK_1) )
    for key in mapK_1:
        if key in uniqMapK:
            assert len( mapK[key] ) == 1
            kmer, c = mapK[key][0]
            if kmer[:-1] == key or kmer[1:] == key: # remove head or tail, they are same
                print (kmer, "remove head or tail same with remove at middle")
                continue
            indelPair.append( (kmer, key, c, mapK_1[key]) )
            #else:
            #    print ("more than one kmer corr. k-1mer", key, mapK_1[key], mapK[key])
    ''' 
    fout = open("indel_pair", "w")
    sortedKmers = sorted(indelKmers)
    for (k1, k2, c1, c2) in sortedKmers:
        fout.write(">kmer_indel%s_1_cov_%s\n" % (ID, c1))
        fout.write("%s\n" % ( k1 ) )
        fout.write(">kmer_indel%s_2_cov_%s\n" % (ID, c2))
        fout.write("%s\n" % ( k2 ) )
        ID += 1
    fout.close()       
    '''
    return indelPair

def update_map_merge(left_i, left_j, key1, key2, mapMerge):

    k1, cov1 = left_i
    k2, cov2 = left_j
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
    for key in left:
        groupSize = len(left[key])
        for i in range(0, groupSize-1):
            for j in range(i+1, groupSize):
                k1, cov1 = left[key][i]
                k2, cov2 = left[key][j]
                if k1[mid] == k2[mid]:
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
                    update_map_merge(left[key][i], left[key][j], key1, key2, mapMerge)
    return mapMerge

def merge_pair(mapMerge):
    pairSet, nonPair = set(), []  
    for (key1, key2) in mapMerge:
        mlen = len( mapMerge [ (key1, key2) ] )
        for i in range(0, mlen-1):
            for j in range(i+1,mlen):
                Left1, Left2, covL1, covL2 = mapMerge[ (key1, key2) ][i]
                Right1, Right2, covR1, covR2 = mapMerge[ (key1, key2) ][j]
                l = len(key1)
                if Left1.find(key1) == 0 and Right1.find(key1) != 0:
                    assert Right1[-l:] == Left1[:l] and Right2[-l:] == Left2[:l]
                    merge1 = Right1 + Left1[l:]
                    merge2 = Right2 + Left2[l:]
                elif Right1.find(key1) == 0 and Left1.find(key1) != 0:
                    assert Left1[-l:] == Right1[:l] and Left2[-l:] == Right2[:l]
                    merge1 = Left1 + Right1[l:]
                    merge2 = Left2 + Right2[l:]
                else:
                    #print ("one side")
                    continue
                small1, small2 = tools.get_smaller_pair_kmer(merge1, merge2)
                if (small1, small2) not in pairSet:
                    pairSet.add( (small1, small2) )
                    nonPair.append( (small1, small2, covL1, covL2, covR1, covR2) )
    return nonPair        


def find_non_pair_kmer(uniqKmer, k):

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
    mapMerge = build_map_merge(left, k)
    nonPair = merge_pair(mapMerge) 

    '''
    fout = open("non_pair", "w")
    sortedKmers = sorted(nonPair)
    ID = 1
    for (k1, k2, c1, c2, c3, c4) in sortedKmers:      
        fout.write(">kmer_non%s_1_cov_%s_cov_%s\n" % (ID, c1, c3))
        fout.write("%s\n" % ( k1 ) )
        fout.write(">kmer_non%s_2_cov_%s_cov_%s\n" % (ID, c2, c4))
        fout.write("%s\n" % ( k2 ) )
        ID += 1
    fout.close()        
    '''
    return nonPair


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
    

def extend_pair(hetePairs, left_index, right_index, k):
    extendPair, extendSet = [], set()
    for (h1, h2, ID) in hetePairs:
        temp1, add1 = extend_to_right(h1, left_index, right_index, k)
        temp2, add2 = extend_to_right(h2, left_index, right_index, k)
        minl = min (len(add1), len(add2) )
        if add1[:minl] != add2[:minl]:
            print (add1, add2)
        #if hamming_distance
        ekmer1, add1 = extend_to_right(temp1, left_index, right_index, k)
        ekmer2, add2 = extend_to_right(temp2, left_index, right_index, k)
        
        minl = min (len(add1), len(add2) )
        if add1[:minl] != add2[:minl]:
            print (add1, add2)
        small1, small2 = tools.get_smaller_pair_kmer(ekmer1, ekmer2)
        if (small1, small2) not in extendSet:
            extendSet.add( (small1, small2) )
            extendPair.append( (small1, small2, ID) )
    return extendPair   



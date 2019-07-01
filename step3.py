#########################################################################
# File Name: step1.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Wed 08 May 2019 11:05:49 AEST
#########################################################################
#!/bin/bash
import sys
import tools
# shift k-mer


def read_unique_kmer(unique_filename):

    left_unique_index = {}
    right_unique_index = {}
    count = 0
    with open(unique_filename, "r") as f:
        for line in f:
            count += 1
            words = line.strip().split()
            key = words[0][:-1]
            if key not in left_unique_index:
                left_unique_index[key] = list()
            left_unique_index[key].append(words[0][-1])

            key = words[0][1:]
            if key not in right_unique_index:
                right_unique_index[key] = list()
            right_unique_index[key].append(words[0][0])
            if count%100000 == 0:
                print "read unique kmer: ", count


    return left_unique_index, right_unique_index

def shift_kmer(h, left_unique_index, right_unique_index, label): 
    # last para forward string or backward string
   
    group = []
    #Rgroup = []

    temp = h
    length = len(h)/2
    #########################
    #x:  AAAATC....TTATT
    #y:   AAATC....TTATT
    # x<x' but y>y'
    #x':  AATA ...TTTT 
    #y': AAATA ...TTT
    #if unique kmer only store one samller kmer this code wrong
    #unique kmer now store forward and bacward string to use this code
    #we later improve this
    #now code is suit for only store one smaller kmer
    ###################

    if label == "l":
        key = temp[1:]
        Rkey = tools.reverse(key)
        #print "ls", h, key, Rkey
        for i in range(0, length):
            flag, flagR = False, False
            if key in left_unique_index and len(left_unique_index[key]) == 1:
                temp = key + left_unique_index[key][0]
                flag = True

            if Rkey in right_unique_index and len(right_unique_index[Rkey]) == 1:
                Rtemp = right_unique_index[Rkey][0] + Rkey
                flagR = True

            if flag == True and flagR == True:

                #print "shift stop, size of group", len(group)
                if len(group) >= 5:
                    print "shift stop", key, temp, "orginal", h 
                    print "forward and backward both have next"
                    print Rkey, Rtemp
                break
                #sys.exit()
            elif flag == True and flagR == False:
                key = temp[1:]
                Rkey = tools.reverse(key)
                group.append( (temp, 'f') )
                #print i, temp, key, Rkey
            elif flag == False and flagR == True:  
                Rkey = Rtemp[:-1]
                key = tools.reverse(Rkey)
                #ward = tool.reverse_ward(ward)
                group.append( (Rtemp, 'b') )
                #print i, Rtemp, key, Rkey
            else:

                #print "shift stop, size of group", len(group)
                if len(group) >= 5:
                    print "shift stop", key, Rkey, "orginal", h
                    print key in left_unique_index, Rkey in right_unique_index
                break
        #print "le", group

    if label == "r":
        key = temp[:-1]
        Rkey = tools.reverse(key)
        #print "rs", h, key, Rkey
        for i in range(0, length):
            flag, flagR = False, False

            if key in right_unique_index and len(right_unique_index[key]) == 1:
                temp = right_unique_index[key][0] + key
                flag = True

            if Rkey in left_unique_index and len(left_unique_index[Rkey]) == 1:
                Rtemp = Rkey + left_unique_index[Rkey][0]
                flagR = True

            if flag == True and flagR == True:
                #print "rr", key, temp
                #print Rkey, Rtemp
                break
            elif flag == True and flagR == False:
                key = temp[:-1]
                Rkey = tools.reverse(key)
                group.append( (temp, 'f') )
                #print i, temp, key, Rkey
            elif flag == False and flagR == True:  
                Rkey = Rtemp[1:]
                key = tools.reverse(Rkey)
                #ward = tool.reverse_ward(ward)
                group.append( (Rtemp, 'b') )
                #print i, Rtemp, key, Rkey
            else:
                break
        #print "re", group    

    ''' 
    if label == "r":
        for i in range(0, length):
            key = temp[:-1]
            if key in unique_index and len(right_unique_index[key]) == 1:
                temp = right_unique_index[key][0] + key
                group.append( temp )
            else:
                break
    '''            
    assert len(group) <= len(h)
    return group 

def group_shift_kmer(pair_filename, unique_filename, output_filename, NGS_kmer):

    left_unique_index, right_unique_index = read_unique_kmer(unique_filename)

    fout = open(output_filename, "w")
    count = 0
    kmers={}
    with open(pair_filename, "r") as f:
        for line in f:
            words = line.strip().split()
            h1 = words[0]
            h2 = words[1]
            leftGroup1 = shift_kmer(h1, left_unique_index, right_unique_index, "l")
            leftGroup2 = shift_kmer(h2, left_unique_index, right_unique_index, "l")
           
            cov1 = int(words[2])
            cov2 = int(words[3])
            sum_coverage = cov1 + cov2 
                
            Flag = True
            i=0
            lenG = min ( len(leftGroup1), len(leftGroup2) )
            while i < lenG:
                if ( tools.hamming_distance(leftGroup1[i][0], leftGroup2[i][0]) != 1 
                   and tools.hamming_distance(leftGroup1[i][0], tools.reverse(leftGroup2[i][0])) != 1 ):
                    Flag = False
                    break
                i += 1
                
            '''
            # 60 equal to lowest coverage
            if cov1 - cov2 >= 60 or cov2 - cov1 >= 60: # or float(cov2/cov1) >=1.7 or float(cov1/cov2) >=1.7:
                continue
            if sum_coverage <= 180 or sum_coverage >= 280:
                continue
            '''    
            if Flag == False:
                continue
               
            rightGroup1 = shift_kmer(h1, left_unique_index, right_unique_index, "r")
            rightGroup2 = shift_kmer(h2, left_unique_index, right_unique_index, "r")
             
            i=0
            lenG = min ( len(rightGroup1), len(rightGroup2) )
            while i < lenG:
                if ( tools.hamming_distance(rightGroup1[i][0], rightGroup2[i][0]) != 1 
                   and tools.hamming_distance(rightGroup1[i][0], tools.reverse(rightGroup2[i][0])) != 1 ):
                    Flag = False
                    break
                i += 1
            if Flag == False:
                continue
            
            group1, group2 =[ (h1, 'f') ], [ (h2, 'f') ]
            group1.extend(leftGroup1) 
            group2.extend(leftGroup2)
            group1.extend(rightGroup1) 
            group2.extend(rightGroup2)

            gSize1 = len(group1)
            gSize2 = len(group2)
            if gSize1 <=kmerSize-1 or gSize2 <= kmerSize-1:
            #if gSize1 <=kmerSize/2 or gSize2 <= kmerSize/2: 
                continue
            count += 1
            fout.write("group %s %s %s %s %s\n" % (count, h1, h2, words[2], words[3]) )
            print "group", count, len(leftGroup1), len(rightGroup1), len(leftGroup2), len(rightGroup2)
            for ele in group1: 
                fout.write( "%s %s " % (ele[0], ele[1]) )
                if ele[0] not in kmers:
                    kmers[ ele[0] ] = []
                kmers[ ele[0] ].append(str(count)+ele[1]+'A') # A is zore   
            fout.write("\n")

            for ele in group2: 
                fout.write("%s %s " % (ele[0] , ele[1] ) )
                if ele[0] not in kmers:
                    kmers[ ele[0] ] = []
                kmers[ ele[0] ].append(str(count) + ele[1] + 'B') # B is one
            fout.write("\n")
    fout.close()
    fout = open(NGS_kmer, "w")
    print "total group number", count
    sortedKmers = sorted(kmers.items())
    filterGroup = set()
    for ele in sortedKmers:
        if len(ele[1]) >= 2:
            print ele
            for ID in ele[1]:
                filterGroup.add(ID[:-2])
            continue
        fout.write("%s" % ele[0])
        l = ele[1]
        for e in l: 
            fout.write(" %s" % e) 
        fout.write("\n")
    fout.close()
    
    print "filter group size", len(filterGroup)
    foutFilter = open(filter_filename, "w") 
    with open(output_filename, "r") as f:
        state = 0
        for line in f:
            if state == 0 and line.startswith("group"):
                words = line.split()
                if words[1] not in filterGroup:
                    foutFilter.write(line)
                    state = 1
                else:
                    state = -1
            elif state == 1:
                foutFilter.write(line)
                state = 2
            elif state == 2:
                foutFilter.write(line)
                state = 0
            elif state == -1:
                state = 0
        


unique_filename=  sys.argv[1] + ".unique.kmer"
pair_filename=  sys.argv[1] + ".pair.kmer"
output_filename= sys.argv[1] + ".group.kmer"
filter_filename= sys.argv[1] + ".filter.group.kmer"
NGS_kmer= sys.argv[1] + ".NGS.kmer"  # this output is sorted
kmerSize = int(sys.argv[2])
group_shift_kmer(pair_filename, unique_filename, output_filename, NGS_kmer)

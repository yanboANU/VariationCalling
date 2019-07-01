#########################################################################
# File Name: compareRealRair.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Thu 09 May 2019 11:51:38 AEST
#########################################################################
#!/bin/bash
import sys
import tools


def read_real_groupID(filename):
    groupID = set()
    with open(filename, "r") as f:
        for line in f:
            words = line.strip().split()
            groupID.add( int(words[0]) ) 
    return groupID 

def read_group_average_weight(filename):
    
    weight = {}
    with open(filename, "r") as f:
        for line in f:
            words= line.strip().split()
            weight[ int(words[0]) ] = round(float(words[-1]), 4)
            '''
            w = round(float(words[-1]), 4)
            if w > 1.01:
                weight[ int(words[0]) ] = w
            '''    
    return weight 


realGroupID = read_real_groupID("GroupID2RefPos")

#weight = read_group_average_weight("deal_matrix.log")
#weight = read_group_average_weight("deal_matrix.log3") 
#weight = read_group_average_weight("groupID_connect_num")

weight = read_group_average_weight("deal_matrix.log4") 
#wkey = set(weight.keys())
#print len(realGroupID), len(wkey), len(wkey.intersection(realGroupID))


#foutTP = open("TPRelateColumn", "w")
#foutFP = open("FPRelateColumn", "w")

foutTP = open("TPAW3_2", "w")
foutFP = open("FPAW3_2", "w")
for ID in weight:
    if ID in realGroupID:
        foutTP.write("%s\n" % weight[ID])
    else:
        foutFP.write("%s\n" % weight[ID])
foutTP.close()
foutFP.close()


#preGroupID = read_real_groupID("aa")
#print len( preGroupID.intersection(realGroupID) )

#!/usr/bin/env/ python3

import sys
import os
import pysam
import string
import copy
import time
import contig
import tools


           #   7    8
best_slope = [1.4, 1.4, 1.7, 1.7, 2.1, 2.1, 2, 2.4, 2.4, 2.7, 2.7, 3, 2.8, 2.8, 2.8, 2.8, 3.1, \
3, 3, 3.1, 3.1, 3.3, 3.3, 3.2, 3.1, 3.3, 3.4, 3.3, 3.4]

class Column:
    # nucleotide: in reference
    # mapContent is a map 'A': support read list
    
    def __init__(self, refPos, a1=0.7, a2=0.3):
        self._ref_pos = refPos

        self._map_content = {}
        self._nucleotide = ''
        self._cov = 0
        self._diff_cov = []
        self._insert_content = {}
        self._insert_number = 0
        self._is_insert = -1
        self._is_mutation = -1
        self._is_delete = -1
        self._is_stable = -1    
        self._a1 = a1
        self._a2 = a2
   
    def _print(self):
        print (self._ref_pos, self._cov)
        print (self._diff_cov)

    def _is_Insert(self):
        
        insertNum = 0
        insertLength = 0
        
        countATCG ={}
        countATCG['A'] = 0
        countATCG['T'] = 0
        countATCG['C'] = 0
        countATCG['G'] = 0
        four = countATCG.keys() 
        
        #if self._cov <7:
            #return False
        for (content, reads) in self._insert_content.items():
            lenReads = len(reads)
            insertNum += lenReads
            insertLength += len(content)*lenReads
            for ele in four:
                if content.count(ele) > 0:
                    countATCG[ele] += lenReads
        '''
        s = tools.sorted_Map_Value_Len(self._insert_content)
        self._insert_content = s   
        # waiting statics data 
        if insertNum > self._cov*0.2 or insertLength > self._cov*0.5:
            #print ("insert")
            return True
        '''
        items = countATCG.items()
        maxb = 0
        maxa = '' 
        for (a,b) in items:
            if b > maxb:
                maxb=b
                maxa=a
        #print (maxb, self._cov, self._get_threshold_for_insert())
        
        if maxb >= self._get_threshold_for_insert():
        #if maxb >= 0.5:
            self._insert_content = maxa
            self._insert_number = maxb
            return True
    
        '''
        if maxb >= self._cov*self._a2:
            self._insert_content = maxa
            self._insert_number = maxb
            return True
        '''    
        return False
    
    def _get_threshold_for_insert(self):
        if self._cov <= 8:
            return self._cov*(0.8)
        if self._cov <= 13:  
            return self._cov*(0.69)
        if self._cov <= 17:  
            return self._cov*(0.55)
        #return 0.5
        
        if self._cov <= 22:
            return self._cov*(0.43)
        if self._cov <= 27:
            return self._cov*(0.39) 
        if self._cov <= 33:
            return self._cov*(0.35) 
        if self._cov <= 37:
            return self._cov*(0.32) 
        if self._cov <= 43:
            return self._cov*(0.3) 
        if self._cov <= 53:
            return self._cov*(0.29) 
        return self._cov*(0.28)
    

    def _get_threshold_for_delete(self):
        if self._cov <= 8:
            return self._cov*(0.7)
        if self._cov <= 12:  
            return self._cov*(0.56)
        #return 0.5
        
        if self._cov <= 17:  
            return self._cov*(0.47)
        if self._cov <= 22:
            return self._cov*(0.38)
        if self._cov <= 27:  
            return self._cov*(0.33)
        if self._cov <= 32:
            return self._cov*(0.29) 
        if self._cov <= 37:  
            return self._cov*(0.28)
        if self._cov <= 42:
            return self._cov*(0.27) 
        if self._cov <= 47:
            return self._cov*(0.25) 
        if self._cov <= 52:
            return self._cov*(0.24) 
        return self._cov*(0.23)
        

    #mutation
    #stable
    #delete

    def _sort_map_content(self): #by value length 
        s = tools.sorted_Map_Value_Len(self._map_content)
        self._map_content = s  
        return s

    def _set_Mutation_or_Delete(self):
        
        self._is_delete = 0
        self._is_mutation = 0
        s = self._sort_map_content()

        #print (self._ref_pos)
        #print (self._map_content)
        #print (s) 
        ''' 
        if self._cov <= 10:
            return 
        '''    
        # do not use it

        #if len(s[0][1]) >= self._cov*self._a1 and s[0][0] != self._nucleotide:
            #print ("reference position %s wrong" % (self._ref_pos))
            #sys.exit("0.7 cov different from ref")
    
        # the latest way to find SNPs
        #print (len(best_slope))
        #print (35-7)
        #sys.exit()
        
        lenS= len(s)

        if lenS ==1 and s[0][0] == '*':
            self._is_delete = 1
        elif lenS == 1:
            return

        if lenS >=2:
            if self._cov >=7:
                if len(s[0][1]) <= best_slope[min(self._cov,35)-7]*len(s[1][1]):
                    if s[0][0] != '*' and s[1][0] != '*':
                        self._is_mutation = 1
                        return

        #for i in range(lenS):
        for i in range(0,2): 
            if s[i][0] == '*' and len(s[i][1]) >= self._get_threshold_for_delete():
            #if s[i][0] == '*' and len(s[i][1]) >= 0.5:
                self._is_delete = 1
                break
            '''
            if s[i][0] == '*' and len(s[i][1]) >= self._cov*self._a1:
                self._is_delete = 1
                break
            ''' 
        ''' 
        # first round of guofei 
        if len(s) >=2 and len(s[0][1]) <= 2*len(s[1][1]):
            if s[0][0] == '*' or s[1][0] == '*':
                self._is_delete = 1
            else:
                self._is_mutation = 1 
        
        
        # second round of guofei
        if len(s)>=2 and len(s[0][1]) <= 3*len(s[1][1]):
            if s[0][0] == '*' or s[1][0] == '*':
                self._is_delete = 1
            else:
                self._is_mutation = 1 
        
        # our original method 
        if len(s[0][1]) < self._cov*self._a1 and len(s[1][1]) >= self._cov*self._a2:
            if s[0][0] == '*' or s[1][0] == '*':
                self._is_delete = 1
            else:
                self._is_mutation = 1 
        '''      
 
    def _set_Lable(self):
        
       
        self._set_Mutation_or_Delete()
      
        if self._is_Insert():
            #print ("insert")
            self._is_insert = 1  
        else:
            self._is_insert = 0 

        if self._is_insert ==0 and self._is_mutation == 0 and self._is_delete == 0:
            self._is_stable = 1
        else:
            self._is_stable = 0 
        #print (self._is_insert, self._is_mutation, self._is_delete, self._is_stable)    

def init_Columns(bamfile, contig, writeLable=False, start=0, end=1000000, a1=0.3, a2=0.4):
    
    columns = {}
    time1 = time.clock()
    '''  
    for read in bamfile.fetch(contig._name, start, end):
        print (read.query_name)
        print (read.reference_name)
        sys.exit()
    '''
    print ("chr name", contig._name)

    ################################################################
    # position start with 0 or 1 ?
    # in the alignedPairs, reference start with 0, read start with 0
    # samtool tview first line, reference position label start with 1 
    ################################################################### 
    #for read in bamfile.fetch(contig._name, start, end):
    #for read in bamfile.fetch("1", start, end):

    for read in bamfile.fetch():
        alignedPairs = read.get_aligned_pairs()
        readName = read.query_name

        print (readName)
        sys.exit()
        i = 0
        label = 1
        lenAlignedPairs = len(alignedPairs) 
        while i+1 < lenAlignedPairs:
            (readPos, referPos) = alignedPairs[i] #check referPos always increse or not
            # "get reference first not None"
            print (i, readPos, referPos)
            #sys.exit()
            while label and isinstance(readPos, int) and not type(referPos) == int:
                i += 1
                (readPos, referPos) = alignedPairs[i]  
            label = 0
            # get columns(align info)
            if type(referPos) == int:
                if referPos not in columns:
                    columns[referPos] = Column(referPos, a1, a2)
                columns[referPos]._cov += 1
                columns[referPos]._nucleotide = contig._seq[referPos] # in the contig seq, rreference start with 0
                if type(readPos) == int:
                    c = read.query_sequence[readPos] 
                    if c not in columns[referPos]._map_content:
                        columns[referPos]._map_content[c] = []
                    columns[referPos]._map_content[c].append(readName)
                elif not type(readPos) == int:
                    c = '*'
                    if c not in columns[referPos]._map_content:
                        columns[referPos]._map_content[c] = []
                    columns[referPos]._map_content[c].append(readName)
                i += 1
            
            #get insert info
            if isinstance(readPos, int) and not type(referPos) == int:
                (readPos_pre, referPos_pre) = alignedPairs[i-1]    
                assert type(referPos_pre) == int 
                if referPos_pre >= read.reference_end-1:
                    break 
                #columns[referPos_pre]._insert_content = [] 
                insertSeq = read.query_sequence[readPos] 
                #while i+1 < len(alignedPairs):
                while i+1 < lenAlignedPairs: 
                    i = i+1
                    (readPos, referPos) = alignedPairs[i]
                    if isinstance(readPos, int) and not type(referPos) == int:
                        insertSeq = insertSeq + read.query_sequence[readPos] 
                    else:
                        break
                if referPos_pre not in columns:
                    columns[referPos_pre] = Column(referPos_pre)
                #print ("pos: %s, content: %s" % (referPos_pre, insertSeq))
                if insertSeq not in columns[referPos_pre]._insert_content:
                    columns[referPos_pre]._insert_content[insertSeq] = []

                columns[referPos_pre]._insert_content[insertSeq].append(readName)
        #print ("check", columns[13]._insert_content)        
        #sys.exit()    
        #i = i+1    
    time2 = time.clock()
    
    print ( "init columns(two loops) running %s Seconds" % (time2 - time1) )
    #for speed, large case not write
    if writeLable:
        write_Columns(contig._name+ "_" + str(start)+"_"+str(end) , columns, end-start)
        write_Ins(contig._name+ "_" + str(start)+"_"+str(end) , columns, end-start)
    write_Cov(contig._name+ "_" + str(start)+"_"+str(end) , columns, end-start)        
    return columns


def write_Cov(contigName, columns, length):

    # in the write, reference start with 1
    filename = contigName + "_cov"
    fout = open(filename, "w")
    fout.write("pos coverage nucleotide\n")
    for (refPos, column) in columns.items():
        if refPos <10000 or refPos>length-10000:
            continue
        fout.write("%s %s %s\n" % (refPos+1, column._cov, column._nucleotide))
    fout.close()

def write_Columns(contigName, columns, length):
    # in the write, reference start with 1
    filename = contigName + "_columns"
    fout1 = open(filename, "w")
    for (refPos, column) in columns.items():
        if refPos <10000 or refPos>length-10000:
            continue
        fout1.write("ref pos: %s  cov: %s nucleotide: %s \n" % (refPos+1, column._cov, column._nucleotide))
        tools.write_Map_Len(fout1, column._map_content)
    fout1.close()
   
def write_Ins(contigName, columns, length):

    # in the write, reference start with 1
    filename = contigName + "_Ins"
    fout2 = open(filename, "w")
    for (refPos, column) in columns.items():
        if len(column._insert_content) != 0:
            if refPos <10000 or refPos>length-10000:
                continue
            fout2.write("reference position: %s  coverage: %s reference nucleotide: %s \n" % (refPos+1, column._cov, column._nucleotide))
            tools.write_Map_Len(fout2, column._insert_content)
    fout2.close()



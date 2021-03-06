#########################################################################
# File Name: ../step4_mutiple_thread.py
# Author: yanbo
# mail: liyanbotop@163.com
# Created Time: Tue 14 May 2019 16:36:45 AEST
#########################################################################
#!/bin/bash
import threading
import sys
import Queue
import tools
from multiprocessing import Process
# https://blog.csdn.net/zx8167107/article/details/81083249

def binarySearch(array, t):
    low = 0
    height = len(array)-1
    while low < height:
        mid = (low+height)/2
        if array[mid][0] < t:
            low = mid + 1
        elif array[mid][0] > t:
            height = mid - 1
        else:
            return array[mid]

    return -1

def find_intersection_kmer(TGS, NGS):
    
    intersection = []
    iN, iT = 0,0
    lenT, lenN = len(TGS), len(NGS)

    while iN < lenN and iT < lenT:
        if (NGS[iN][0] < TGS[iT]):
            iN +=1
        elif NGS[iN][0] > TGS[iT]:   
            iT +=1
        elif NGS[iN][0] == TGS[iT]:
            intersection.append( NGS[iN] )
            iN +=1
            iT +=1
        else:
            print ("error")
            iN +=1
            iT +=1
    return intersection       

def decide_kmer(intersection):
    #kmer direction and kmer group A or group B
    #each read for one kmer  one direction and one group
    label = {}
    order = {}

    #seqOrder = 0
    for (ele, pos) in intersection:
        key = ele[:-2]
        l = ele[-2:]

        if key not in order:
            order[ key ] = pos
            #seqOrder += 1 
        if key not in label:
            label[ key ] = {}
            if l not in label[key]:
                label[key][ l ] = 0
            label[key][ l ] += 1
 
    Pos = {}

    for key in label:

        #assert len((label[key])) == 1
        maxSupport = 0
        maxLabel = ""
        values = sorted( label[key].values() )
        if values.count(values[-1]) >= 2:
            # more than one most support, consider a coincidence
            print intersection
            continue

        # only keep kmer have one most support
        for l in label[key]:
            if label[key][l] > maxSupport:
                maxLabel = l
                maxSupport = label[key][l]
        
        if maxLabel[-1] == 'A':       
            Pos[ key ]  = (maxLabel[0] , 0)

        if maxLabel[-1] == 'B':       
            Pos[ key ]  = (maxLabel[0] , 1)


    PosList = []

    orderList = sorted(order.items(), key=lambda item:item[1])
    for (k, o) in orderList: # according kmer happen in TGS reads order
        PosList.append( (k+Pos[k][0], Pos[k][1], o) )
    #print PosList
    return PosList       
    


#class Reader(threading.Thread):
class Reader(Process):
    def __init__(self, thread_id, NGS, temp_list):
        super(Reader, self).__init__()
        self.thread_id = thread_id
        self.NGS = NGS
        self.temp_list = temp_list
 
    def run(self):
        state = 0
        fout = open("part_matrix_"+ str(self.thread_id), "w")
        count = 0
        for text in self.temp_list:
            #columns = text.split('\001')
            if state==0 and text.startswith("@"):
                readID=text.strip()[1:]
                state = 1
                count += 1
                if count % 1000 == 0:
                    print ("thread ", self.thread_id, "deal reads ", count) 
            elif state == 1:
                #print readID
                state = 0
                seq = text.strip() 
                seqLen = len(seq)
                intersection = []
                for i in range(seqLen-21):
                    key = str(seq[i:i+21])
                    Rkey = tools.reverse(key)
                    Rflag = False 
                    if key > Rkey:
                        key = Rkey
                        Rflag = True
                    re = binarySearch(NGS, key)
                    if re != -1:
                        if Rflag == True:
                            #print "1"
                            temp = re[1][:-2] + tools.reverse_ward(re[1][-2]) + re[1][-1] # Rkmer direction is opposite with reads
                            intersection.append( (temp, i) ) # i is position in reads
                        else:    
                            #print "2"
                            intersection.append( (re[1], i) ) # i is position in reads
                            
                #print intersection
                '''
                Rseq = tools.reverse(seq)
                intersection = []
                for i in range(seqLen-21):
                    key = str(Rseq[i:i+21])
                    re = binarySearch(NGS, key)
                    if re != -1:
                        intersection.append( re[1] )
                print intersection
                '''

                if len(intersection) > 0:
                    print intersection
                PosList = decide_kmer(intersection)
                #if count == 10:
                    #sys.exit()
                if len(PosList) <= 1:
                    continue
                fout.write( "%s %s " % ( len(PosList), readID ) )
                for (p, binary, pos) in PosList:
                    fout.write( "%s %s %s " % (p, binary, pos) )
                score=len(PosList)*'4'    
                fout.write("%s\n" % score)
                #print len(seq), len(kmers)
                #sys.exit()
            else:
                continue
        fout.close()
 
 
class Partition(object):
    def __init__(self, file_name, thread_num):
        self.file_name = file_name
        self.block_num = thread_num
 
    def part_and_queue(self):
        pos_list = []
        file_size = tools.file_lines(self.file_name)
        block_size = file_size / self.block_num
        start_pos = 0
        global q
 
        for i in range(self.block_num):
            if i == self.block_num - 1:
                end_pos = file_size - 1
                pos_list.append((start_pos, end_pos))
                break
            end_pos = start_pos + block_size - 1
            if end_pos >= file_size:
                end_pos = file_size - 1
            if start_pos >= file_size:
                break
            pos_list.append((start_pos, end_pos))
            start_pos = end_pos + 1
 
        fd = open(self.file_name, 'r')
        for pos_tu in pos_list:
            temp_text = []
            start = pos_tu[0]
            end = pos_tu[1]
            while start <= end:
                text = fd.readline().strip('\n')
                temp_text.append(text)
                start = start + 1
 
            q.put(temp_text)
        fd.close()
 
 
if __name__ == '__main__':
    
    #TGS_filename = "/media/yanbo/Data/NA12878/46xPacbio_alignment/chr22.fastq"
    #NGS_filename="chr22.NGS.kmer"
    if len(sys.argv) < 3:
        print ("chr22.NGS.kmer, ~/bio/HA/HapCUT2/pacbio/hg19/flag0flag16/chr22.fastq, thred_num")
        sys.exit()

    NGS_filename=sys.argv[1]
    TGS_filename = sys.argv[2]
    
    NGS = []
    with open(NGS_filename, "r") as f:
        for line in f:
            NGS.append(line.strip().split())


    #plat_form = sys.argv[2]
    #a_o_b = sys.argv[3]
 
    #thread number
    thread_num = int(sys.argv[3])
    q = Queue.Queue()
    p = Partition(TGS_filename, thread_num)
    t = []
    p.part_and_queue()

    #global q
    for i in range(thread_num):
        temp_list = q.get()
        t.append(Reader(i, NGS, temp_list))
    for i in range(thread_num):
        t[i].start()
    for i in range(thread_num):
        t[i].join()



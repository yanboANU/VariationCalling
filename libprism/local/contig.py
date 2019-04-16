#!/usr/bin/env python3


class Contig:
    def __init__(self, name, seq):
        self._name = name
        self._seq = seq
        self._len = len(seq) 
     
    def print(self):
        print (self._name, self._len)

def read_Contig(filename):
    contigs = {}
    f = open(filename, "r")
    
    seq = ""
    for line in f:
        if line.startswith(">"):
            if len(seq) > 0:
                contig = Contig(name,seq)  
                contig.print()
                contigs[name] = contig
            name = line.split()[0][1:]
            seq = "" 
        else:
            seq = seq + line.strip()
            
    if len(seq) > 0:
        contig = Contig(name,seq)  
        contig.print()
        contigs[name] = contig
    f.close()        
    return contigs

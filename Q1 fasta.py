import matplotlib.pyplot as plt
import numpy
import re

#Useful Methods

#Parse through fasta data to separate details from sequence.
#Returns tuple of name and sequence

def parse_fasta(fasta):
    seq = ""
    f = open(fasta, 'r')
    next(f)
    for line in f:
        seq += line.rstrip()
    return seq

#Create a list containing genes.
def GetGenes(GTuples, sequence):
    geneslist = list()
    for GeneTuple in GTuples:
        print GeneTuple
        geneslist.append(sequence[int(GeneTuple[0]):int(GeneTuple[1])])
    return geneslist

#Make list of ratio of GC bases.
def GCContent(geneslist,lengthslist):
    GCList = []
    for gene, length in zip(geneslist, lengthslist):
        GC = 0.0
        for base in gene:
            if base == 'G' or base == 'C':
                GC += 1
        GCList.append(GC / length)
    return GCList

#reverse compliment
def reverseC(seq):
    RC = ""
    for base in range(len(seq)-1,-1,-1):
        if seq[base] == "A":
            RC += "T"
        elif seq[base] == "T":
            RC += "A"
        elif seq[base] == "C":
            RC += "G"
        else:
            RC += "C"
    return RC

"""
######################################################################
#Prepare variables for graph
genelengths = list()
genetuples = list()
genes = list()
GC = list()

#Fastadata variable to access for other functions down the line.
fastadata = parse_fasta('s_meliloti.fa')

#Read file by line and populate genelengths and genetuples from GFF3

fh = open('s_meliloti.gff3', 'r')
for line in fh:
    if re.search('gene', line) is not None:
        values = line.split('\t')
        start = int(values[3]) #start coordinate
        end = int(values[4]) #end coordinate
        genelength = end - start + 1
        genetuples.append((start,end))
        genelengths.append(genelength)

genes = GetGenes(genetuples, fastadata)
GC = GCContent(genes, genelengths)
plt.scatter(genelengths, GC)
plt.show()
"""

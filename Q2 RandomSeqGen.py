import random
import re

def seqCompare(seq1, seq2, detailed = False):
    seqCompare = ''
    for letter1, letter2 in zip(seq1, seq2):
        if letter1 != letter2:
            seqCompare += '-'
        else:
            seqCompare += letter1
    if detailed:  
        print('Character comparison')
        print('Dashes indicate different positions. Nucleotides seen are conserved between sequences.')
        print('-'*50)
        print('Given')
        print(seq1)
        print('-'*50)
        print('Conserved')
        print(seqCompare)
        print('-'*50)
        print('Generated')
        print(seq2)
    return(seqCompare)

def randomSeqGen(seqLength, pA, pT, pC, pG):
    seqRanGen = ''
    count = 0
    while count < seqLength:
        RNG = random.random()
        if RNG <= pA:
            seqRanGen += 'A'
        elif RNG > pA and RNG <= pA+pC:
            seqRanGen += 'C'
        elif RNG > pA+pC and RNG <= pA+pC+pG:
            seqRanGen += 'G'
        else:
            seqRanGen += 'T'
        count += 1
    return seqRanGen

def multiSeqGen(number, lengths):
    sequences = []
    count = 0
    while count < number:
        sequences.append(randomSeqGen(lengths))
        count += 1
    return sequences

def randomRead(sequence, readLength, reads = 1):
    RNG = random.randint(0,len(sequence) - readLength)
    readsList = []
    count1 = 0
    while count1 < reads:
        read = ''
        count = 0
        while count < readLength:
            read += sequence[RNG+count]
            count += 1
        count1 += 1
        readsList.append(read)
    f = open(str(output), 'w')
    for x in readsList:
        f.write(str(x) + '\n')
    return readsList

###Problem - reads same output every time
def randomReadOutput(sequence, readLength, reads = 1, output = None):
    RNG = random.randint(0,len(sequence) - readLength)
    count1 = 0
    f = open(str(output), 'w')
    while count1 < reads:
        read = ''
        count = 0
        while count < readLength:
            read += sequence[RNG+count]
            count += 1
        count1 += 1
        f.write(read + '\n')
    return

def readfile(inputfile):
    sequenceList = []
    f = open(inputfile, 'r')
    for line in f:
        sequenceList.append(line)
    return sequenceList

#function countKmers
#takes list of sequences (sequenceList) and runs seqCompare on all of them, stores in compareList
#compareList is filtered by length kLength, finding only k-mers of kLength parameter
#kdict stores all k-mers found and their frequencies
#kmode finds the mode of kdict if it exists, otherwise it calculates the mean
def countKmers(sequenceList, kLength = 17):
    kdict = {}
    countRow = 0
    countCol = 1
    coverage = 0
    #populate compareList with seqCompare for each sequence comparison only once
    #each comparison's kmers are appended to compareList
    while countRow < len(sequenceList) - 1:
        while countCol < len(sequenceList):
            seqSplit = filter(lambda x: len(x) == kLength, seqCompare(sequenceList[countRow],sequenceList[countCol]).split('-'))
            for kmer in seqSplit:
                if kmer in kdict:
                    kdict[kmer] += 1
                else:
                    kdict[kmer] = 1
            countCol += 1
        countRow += 1
        print 'Row ' + str(countRow) + ' done.'
        countCol = countRow + 1
    mode = 1
    if kdict.values():
        mode = sorted(kdict.values(), key=kdict.values().count)[-1]
    if mode == 1 and kdict.values():
            coverage = float(sum(kdict.values) / len(kdict.values))
    else:
        coverage = mode
    if coverage != 0:
        print 'Predicted coverage with k-length ' + str(kLength) + ' : ' + str(coverage)
        print 'Predicted sequence length: ' + str(len(sequenceList)*len(sequenceList[0])/coverage)
    else:
        print 'Could not find good estimate of coverage'
    return coverage

#Function that gives mutated version of a given sequence with
#  similar distribution of nucleotides
#Get sample sequence and nuc distribution

#single nuc
def randomSeqGenMk1(sequence):
    numA = 0.0
    numT = 0.0
    numC = 0.0
    numG = 0.0
    for letter in sequence:
        if letter == 'A':
            numA +=1
        elif letter == 'T':
            numT += 1
        elif letter == 'C':
            numC += 1
        elif letter == 'G':
            numG += 1
        else:
            print ('Incorrect seqeuence letter found: ' + letter)
            break
    #Determine probabilities of each nuc from distributions

    seqLen = len(sequence)
    probA = numA/seqLen
    probT = numT/seqLen
    probC = numC/seqLen
    probG = numG/seqLen


    #Generate Sequence
    seqRanGen = ''
    count = 0
    while count < seqLen:
        RNG = random.random()
        if RNG <= probA:
            seqRanGen += 'A'
        elif RNG <= (probA+probT) and RNG > probA:
            seqRanGen += 'T'
        elif RNG <= (probA+probT+probC) and RNG > (probA+probT):
            seqRanGen += 'C'
        else:
            seqRanGen += 'G'
        count+= 1

    #Do some stats on the generated sequence to compare to original
    numARNG = 0.0
    numTRNG = 0.0
    numCRNG = 0.0
    numGRNG = 0.0

    for letter in seqRanGen:
        if letter == 'A':
            numARNG +=1
        elif letter == 'T':
            numTRNG += 1
        elif letter == 'C':
            numCRNG += 1
        elif letter == 'G':
            numGRNG += 1
        else:
            print ('Incorrect seqeuence letter found: ' + letter)
            break
    print("Random sequence generated:")
    print(seqRanGen)
    print('-'*20)
    print('Distributions ')
    print('initial    /   generated')
    print('-'*20)
    print('A: ' + str(numA) + ' '*10 + str(numARNG))
    print('T: ' + str(numT) + ' '*10 + str(numTRNG))
    print('C: ' + str(numC) + ' '*10 + str(numCRNG))
    print('G: ' + str(numG) + ' '*10 + str(numGRNG))
    print('')

    #Check similarity between orignal sequence and generated
    #Determined by #chars in same position
    count = 0
    similarities = 0.0
    for nuc in sequence:
        if nuc == seqRanGen[count]:
            similarities += 1.0
        count += 1
    percentSim = similarities/float(seqLen)*100

    seqCompare(sequence,seqRanGen)
    print('-'*50)
    print('Analysis of given and generated...')
    print('Characters in similar place: ' + str(similarities))
    print('Percent similarity: %' + str(percentSim))

    return seqRanGen

def randomSeqGenMk2(sequence):
    numA = 0.0
    numT = 0.0
    numC = 0.0
    numG = 0.0
    dinucDict = {'AA':0, 'AT':0, 'AC':0, 'AG':0,
                 'TA':0, 'TT':0, 'TC':0, 'TG':0,
                 'CA':0, 'CT':0, 'CC':0, 'CG':0,
                 'GA':0, 'GT':0, 'GC':0, 'GG':0}
    prevLetter = ''
    for letter in sequence:
        if prevLetter:
            if letter == 'A':
                numA +=1
                key = prevLetter + 'A'
                dinucDict[key] += 1
                prevLetter = 'A'
            elif letter == 'T':
                numT += 1
                key = prevLetter + 'T'
                dinucDict[key] += 1
                prevLetter = 'T'
            elif letter == 'C':
                numC += 1
                key = prevLetter + 'C'
                dinucDict[key] += 1
                prevLetter = 'C'
            elif letter == 'G':
                numG += 1
                key = prevLetter + 'G'
                dinucDict[key] += 1
                prevLetter = 'G'
            else:
                print ('Incorrect seqeuence letter found: ' + letter)
                break
        else:
            if letter == 'A':
                numA +=1
                prevLetter = 'A'
            elif letter == 'T':
                numT += 1
                prevLetter = 'T' 
            elif letter == 'C':
                numC += 1
                prevLetter = 'C'
            elif letter == 'G':
                numG += 1
                prevLetter = 'G'
            else:
                print ('Incorrect seqeuence letter found: ' + letter)
                break
    #For first letter
    seqRanGen = ''
    prevLetter = ''
    seqLen = len(sequence)
    probA = numA/seqLen
    probT = numT/seqLen
    probC = numC/seqLen
    probG = numG/seqLen
    RNG1 = random.random()
    if RNG1 <= probA:
        seqRanGen += 'A'
        prevLetter = 'A'
    elif RNG1 <= (probA+probT) and RNG1 > probA:
        seqRanGen += 'T'
        prevLetter = 'T'
    elif RNG1 <= (probA+probT+probC) and RNG1 > (probA+probT):
        seqRanGen += 'C'
        prevLetter = 'C'
    else:
        seqRanGen += 'G'
        prevLetter = 'G'
    #Rest of sequence
    count = 1
    probAA = dinucDict['AA']/numA
    probAT = dinucDict['AT']/numA
    probAC = dinucDict['AC']/numA
    probAG = dinucDict['AG']/numA
    probTA = dinucDict['TA']/numT
    probTT = dinucDict['TT']/numT
    probTC = dinucDict['TC']/numT
    probTG = dinucDict['TG']/numT
    probCA = dinucDict['CA']/numC
    probCT = dinucDict['CT']/numC
    probCC = dinucDict['CC']/numC
    probCG = dinucDict['CG']/numC
    probGA = dinucDict['GA']/numG
    probGT = dinucDict['GT']/numG
    probGC = dinucDict['GC']/numG
    probGG = dinucDict['GG']/numG
    while count < seqLen:
        RNG1 = random.random()
        if prevLetter == 'A':
            if RNG1 <= probAA:
                seqRanGen += 'A'
                prevLetter = 'A'
            elif RNG1 <= (probAA+probAT) and RNG1 > probAA:
                seqRanGen += 'T'
                prevLetter = 'T'
            elif RNG1 <= (probAA+probAT+probAC) and RNG1 > (probAA+probAT):
                seqRanGen += 'C'
                prevLetter = 'C'
            else:
                seqRanGen += 'G'
                prevLetter = 'G'
        elif prevLetter == 'T':
            if RNG1 <= probTA:
                seqRanGen += 'A'
                prevLetter = 'A'
            elif RNG1 <= (probTA+probTT) and RNG1 > probTA:
                seqRanGen += 'T'
                prevLetter = 'T'
            elif RNG1 <= (probTA+probTT+probTC) and RNG1 > (probTA+probTT):
                seqRanGen += 'C'
                prevLetter = 'C'
            else:
                seqRanGen += 'G'
                prevLetter = 'G'
        elif prevLetter == 'C':
            if RNG1 <= probCA:
                seqRanGen += 'A'
                prevLetter = 'A'
            elif RNG1 <= (probCA+probCT) and RNG1 > probCA:
                seqRanGen += 'T'
                prevLetter = 'T'
            elif RNG1 <= (probCA+probCT+probCC) and RNG1 > (probCA+probCT):
                seqRanGen += 'C'
                prevLetter = 'C'
            else:
                seqRanGen += 'G'
                prevLetter = 'G'
        else:
            if RNG1 <= probGA:
                seqRanGen += 'A'
                prevLetter = 'A'
            elif RNG1 <= (probGA+probGT) and RNG1 > probGA:
                seqRanGen += 'T'
                prevLetter = 'T'
            elif RNG1 <= (probGA+probGT+probGC) and RNG1 > (probGA+probGT):
                seqRanGen += 'C'
                prevLetter = 'C'
            else:
                seqRanGen += 'G'
                prevLetter = 'G'
        count+= 1
    print ''
    print 'Sequence generated randomly with dependence'
    print seqRanGen
    #Compare new to given
    dinucDict2 = {'AA':0, 'AT':0, 'AC':0, 'AG':0,
             'TA':0, 'TT':0, 'TC':0, 'TG':0,
             'CA':0, 'CT':0, 'CC':0, 'CG':0,
             'GA':0, 'GT':0, 'GC':0, 'GG':0}
    for letter in seqRanGen:
        if prevLetter:
            if letter == 'A':
                numA +=1
                key = prevLetter + 'A'
                dinucDict2[key] += 1
                prevLetter = 'A'
            elif letter == 'T':
                numT += 1
                key = prevLetter + 'T'
                dinucDict2[key] += 1
                prevLetter = 'T'
            elif letter == 'C':
                numC += 1
                key = prevLetter + 'C'
                dinucDict2[key] += 1
                prevLetter = 'C'
            elif letter == 'G':
                numG += 1
                key = prevLetter + 'G'
                dinucDict2[key] += 1
                prevLetter = 'G'
            else:
                print ('Incorrect seqeuence letter found: ' + letter)
                break
        else:
            if letter == 'A':
                numA +=1
                prevLetter = 'A'
            elif letter == 'T':
                numT += 1
                prevLetter = 'T' 
            elif letter == 'C':
                numC += 1
                prevLetter = 'C'
            elif letter == 'G':
                numG += 1
                prevLetter = 'G'
            else:
                print ('Incorrect seqeuence letter found: ' + letter)
                break
    print ''
    
    print 'Differences in dinucleotide distribution'
    diff = [k for k in dinucDict if dinucDict[k] != dinucDict2[k]]
    totalDiffs = 0
    for k in diff:
        print k, ':', dinucDict[k], '->', dinucDict2[k]
        totalDiffs += abs(dinucDict[k]-dinucDict2[k])
    print 'Total number of differences: ' + str(totalDiffs)
    print 'Percent similarity: ' + str(100.0 - float(totalDiffs)/(seqLen-1)*100)
    print ''
    seqCompare(sequence, seqRanGen)
    return seqRanGen
    
###main
#a = randomSeqGen(1000000)
#print 'Random sequence generated'
#b = randomRead(a, 100, 10000)
#print 'Random reads generated'
#c = countKmers(b)
maincount = 0
seqList = []
while maincount < 1000:
    seqList.append(randomSeqGen(120,.25,.2,.3,.25))
    maincount += 1
matches = 0
for seq in seqList:
    if 'GCCC' in seq:
        matches += 1.0
print str(matches/float(len(seqList)))

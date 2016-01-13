import re
import math
import random

#reads paramter file and makes probability dictionary
def readPD(filename):
    PD = {}
    with open(filename, 'r') as f:
        for line in f:
            prev = ''
            for word in line.split():
                if re.search('[ACGT]', word):
                    prev = word
                else:
                    PD[prev] = float(word)
    return PD

#turns a dinucleotide count distribution into a probability distrition
def toProb(dinucDict):
    probDict = {'A':0.0,  'T':0.0,  'C':0.0,  'G':0.0,
                'AA':0.0, 'AT':0.0, 'AC':0.0, 'AG':0.0,
                'TA':0.0, 'TT':0.0, 'TC':0.0, 'TG':0.0,
                'CA':0.0, 'CT':0.0, 'CC':0.0, 'CG':0.0,
                'GA':0.0, 'GT':0.0, 'GC':0.0, 'GG':0.0}
    countA = float(dinucDict['AA']) + float(dinucDict['AC']) + float(dinucDict['AG']) + float(dinucDict['AT'])
    countT = float(dinucDict['TA']) + float(dinucDict['TC']) + float(dinucDict['TG']) + float(dinucDict['TT'])
    countC = float(dinucDict['CA']) + float(dinucDict['CC']) + float(dinucDict['CG']) + float(dinucDict['CT'])
    countG = float(dinucDict['GA']) + float(dinucDict['GC']) + float(dinucDict['GG']) + float(dinucDict['GT'])
    for key in probDict:
        if key[0] == 'A' and len(key) > 1:
            probDict[key] = float(dinucDict[key]) / countA
        elif key[0] == 'T' and len(key) > 1:
            probDict[key] = float(dinucDict[key]) / countT
        elif key[0] == 'C' and len(key) > 1:
            probDict[key] = float(dinucDict[key]) / countC
        elif key[0] == 'G' and len(key) > 1:
            probDict[key] = float(dinucDict[key]) / countG
    return probDict

def markovGen(seqLength, PD):
    #probability dictionary
    #PD = {'A':.25, 'C':.25, 'G':.25, 'T':25,
    #            'AA':.4, 'AC':.4, 'AG':.1, 'AT':.1,
    #            'TA':.6, 'TC':.1, 'TG':.1, 'TT':.1,
    #            'CA':.2, 'CC':.2, 'CG':.3, 'CT':.3,
    #            'GA':.25, 'GC':.25, 'GG':.1, 'GT':.4}
    seqRanGen = ''
    prevLetter = ''
    RNG1 = random.random()
    #first letter
    if RNG1 <= PD['A']:
        seqRanGen += 'A'
        prevLetter = 'A'
    elif RNG1 <= (PD['A']+PD['C']) and RNG1 > PD['A']:
        seqRanGen += 'C'
        prevLetter = 'C'
    elif RNG1 <= (PD['A']+PD['C']+PD['G']) and RNG1 > (PD['A']+PD['C']):
        seqRanGen += 'G'
        prevLetter = 'G'
    else:
        seqRanGen += 'T'
        prevLetter = 'T'
    #rest
    while len(seqRanGen) < seqLength:
        RNG1 = random.random()
        if prevLetter == 'A':
            if RNG1 <= PD['AA']:
                seqRanGen += 'A'
                prevLetter = 'A'
            elif RNG1 <= (PD['AA']+PD['AC']) and RNG1 > PD['AA']:
                seqRanGen += 'C'
                prevLetter = 'C'
            elif RNG1 <= (PD['AA']+PD['AC']+PD['AG']) and RNG1 > (PD['AA']+PD['AC']):
                seqRanGen += 'G'
                prevLetter = 'G'
            else:
                seqRanGen += 'T'
                prevLetter = 'T'
        elif prevLetter == 'C':
            if RNG1 <= PD['CA']:
                seqRanGen += 'A'
                prevLetter = 'A'
            elif RNG1 <= (PD['CA']+PD['CC']) and RNG1 > PD['CA']:
                seqRanGen += 'C'
                prevLetter = 'C'
            elif RNG1 <= (PD['CA']+PD['CC']+PD['CT']) and RNG1 > (PD['CA']+PD['CC']):
                seqRanGen += 'G'
                prevLetter = 'G'
            else:
                seqRanGen += 'T'
                prevLetter = 'T'
        elif prevLetter == 'G':
            if RNG1 <= PD['GA']:
                seqRanGen += 'A'
                prevLetter = 'A'
            elif RNG1 <= (PD['GA']+PD['GC']) and RNG1 > PD['GA']:
                seqRanGen += 'C'
                prevLetter = 'C'
            elif RNG1 <= (PD['GA']+PD['GC']+PD['GG']) and RNG1 > (PD['GA']+PD['GC']):
                seqRanGen += 'G'
                prevLetter = 'G'
            else:
                seqRanGen += 'T'
                prevLetter = 'T'
        else:
            if RNG1 <= PD['TA']:
                seqRanGen += 'A'
                prevLetter = 'A'
            elif RNG1 <= (PD['TA']+PD['TC']) and RNG1 > PD['TA']:
                seqRanGen += 'C'
                prevLetter = 'C'
            elif RNG1 <= (PD['TA']+PD['TC']+PD['TG']) and RNG1 > (PD['TA']+PD['TC']):
                seqRanGen += 'G'
                prevLetter = 'G'
            else:
                seqRanGen += 'T'
                prevLetter = 'T'
    return seqRanGen


def markovAnalyze(sequence, modelPP, PDs):
    dinucDict = {'A':.25, 'T':.25, 'C':.25, 'G':.25,
             'AA':0, 'AT':0, 'AC':0, 'AG':0,
             'TA':0, 'TT':0, 'TC':0, 'TG':0,
             'CA':0, 'CT':0, 'CC':0, 'CG':0,
             'GA':0, 'GT':0, 'GC':0, 'GG':0}
    prevLetter = ''
    #populate dinucDict with count for each dinucleotide
    for letter in sequence:
        if prevLetter:
            if letter == 'A':
                key = prevLetter + 'A'
                dinucDict[key] += 1
                prevLetter = 'A'
            elif letter == 'T':
                key = prevLetter + 'T'
                dinucDict[key] += 1
                prevLetter = 'T'
            elif letter == 'C':
                key = prevLetter + 'C'
                dinucDict[key] += 1
                prevLetter = 'C'
            elif letter == 'G':
                key = prevLetter + 'G'
                dinucDict[key] += 1
                prevLetter = 'G'
            else:
                #print ('Incorrect seqeuence letter found: ' + letter)
                break
        else:
            if letter == 'A':
                prevLetter = 'A'
            elif letter == 'T':
                prevLetter = 'T' 
            elif letter == 'C':
                prevLetter = 'C'
            elif letter == 'G':
                prevLetter = 'G'
            else:
                #print ('Incorrect seqeuence letter found: ' + letter)
                break
    #print dinucDict
    probDict = toProb(dinucDict)
    #print probDict
    #compare found PD to given PDs
    difList = []
    difSum = 0.0
    for PD in PDs:
        dif = 0.0
        for key in PD:
            if len(key) != 1:
                dif += math.pow((abs(float(PD[key]) - float(probDict[key]))*100), 2)
        #print dif
        difSum += dif
        difList.append([PD, dif])
    for pair in difList:
        pair[1] = 1 - (pair[1] / difSum)
    #print difList
    #compute posterior probability
    #print 'Posterior probability'
    probs = []
    sumProbs = 0.0
    for PD in PDs:
        prev = ''
        prob = 1.0
        for nuc in sequence:
            if prev:
                dinuc = prev + nuc
                prob = prob * float(PD[nuc])
            prev = nuc
        probs.append([PD,prob])
        sumProbs += float(prob)
        print probs
    #for prob in probs:
        #prob[1] = prob[1]/sumProbs
    return probs

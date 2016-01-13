import random
#import numpy
#import matplotlib.pyplot as plt

def seqGen(length):
    seq = ""
    count = 0
    while count < length:
        RNG = random.random()
        if RNG <= 0.25:
            seq = seq + "A"
        elif RNG <= 0.5:
            seq = seq + "C"
        elif RNG <= 0.75:
            seq = seq + "T"
        else:
            seq = seq + "G"
        count += 1
    return seq

def seqScore(seq, scoreArray):
    #score array [Miss(AT), Hit(CG)]
    count = 0

    #create array of partial sums
    partialSum = []
    prev = 0
    while count < len(seq):
        scoreMod = 0
        #if miss(AT), subtract, else add
        if seq[count] == "A": #or seq[count] == "T":
            scoreMod = scoreArray[0]
        else:
            scoreMod = scoreArray[1]
        prev = prev + scoreMod
        partialSum.append(prev)
        count += 1

    #excursion plot
    peaks = []
    recentLow = 0
    count = 0
    while count < len(partialSum):
        if partialSum[count] < recentLow:
            recentLow = partialSum[count]
            count += 1
        else:
            highScore= 0
            start = count
            while recentLow < partialSum[count] and count < len(partialSum):
                if highScore < partialSum[count] - recentLow:
                    end = count
                    highScore = partialSum[count] - recentLow
                    count += 1
                    if count >= len(partialSum):
                        break
                else:
                    count += 1
                    if count >= len(partialSum):
                        break
            peaks.append([highScore,start,end])
            if count < len(partialSum):
                recentLow = partialSum[count]
            count+= 1
##    f = open('quizOutput.txt', 'w')
##    f.write('Sequence given')
##    f.write(seq)
##    f.write('Partial Sums Array')
##    f.write(partialSum)
##    f.write('Local Peaks')
##    f.write(peaks)
##    f.write('Highest Peak')
##    f.write(max(peaks))
    print 'Sequence given'
    print seq
    print 'Partial Sums Array'
    print partialSum
    print 'Local Peaks'
    print peaks
    print 'Highest Peak'
    print max(peaks)
    #plt.plot(partialSum)
    #plt.show()
    return

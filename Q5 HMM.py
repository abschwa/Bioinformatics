import itertools

#states   E (extracellular), M (membranous)
#output H (hydrophilic)  , T (hydrophobic)
#in general, pEH > pET, and pMH < pMT

def naiveHMM(seq, initE, initM, probEE, probEM, probME, probMM):
    stateSeqs = list(itertools.product(['E','M'], repeat = len(seq)))
    stateProbs = {}
    for ss in stateSeqs:
        prob = 1
        prev = ""
        string = ""
        if ss[0] == "E":
            prob = initE
            prev = "E"
            string = string + "E"
        else:
            prob = initM
            prev = "M"
            string = string + "M"
        count = 1
        while count < len(seq):
            if prev == "E" and ss[count] == "E":
                prob = prob * probEE
                count += 1
                string = string + "E"
            elif prev == "E" and ss[count] == "M":
                prob = prob * probEM
                count += 1
                string = string + "M"
                prev = "M"
            elif prev == "M" and ss[count] == "E":
                prob = prob * probME
                count += 1
                string = string + "E"
                prev = "E"
            else:
                prob = prob * probMM
                count+= 1
                string = string + "M"
        stateProbs[string] = prob
    return stateProbs

#EE is given E, prob E next
#ME is given M, prob E next
#transition matrix (hidden as X * hidden as Y)
#confusion matrix (observed as X * hidden as Y)
def forwardHMM(seq, initE, initM, tMat, cMat):
    #rename transition probs for ease
    EE = tMat[0][0]
    EM = tMat[0][1]
    ME = tMat[1][0]
    MM = tMat[1][1]

    EH = cMat[0][0]
    ET = cMat[0][1]
    MH = cMat[1][0]
    MT = cMat[1][1]
    #initialize arrays
    state1 = [0] * len(seq)
    state2 = [0] * len(seq)
    state1[0] = initE
    state2[0] = initM
    i = 1
    while i < len(seq):
        if seq[i] == 'H':
            state1[i] = (state1[i-1]*EE + state2[i-1]*ME) * EH
            state2[i] = (state1[i-1]*EM + state2[i-1]*MM) * MH
            i += 1
        else:
            state1[i] = (state1[i-1]*EE + state2[i-1]*ME) * ET
            state2[i] = (state1[i-1]*EM + state2[i-1]*MM) * MT
            i+= 1
    print state1
    print state2
    return state1[i-1] + state2[i-1]

def viterbiNum(seq, initE, initM, tMat, cMat):
#rename transition probs for ease
    EE = tMat[0][0]
    EM = tMat[0][1]
    ME = tMat[1][0]
    MM = tMat[1][1]

    EH = cMat[0][0]
    ET = cMat[0][1]
    MH = cMat[1][0]
    MT = cMat[1][1]
    #initialize arrays and string
    state1 = [0] * len(seq)
    state2 = [0] * len(seq)
    state1[0] = initE
    state2[0] = initM
    i = 1
    while i < len(seq):
        if seq[i] == 'H':
            state1[i] = max((state1[i-1]*EE, state2[i-1]*ME)) * EH
            state2[i] = max((state1[i-1]*EM, state2[i-1]*MM)) * MH
            i += 1
        else:
            state1[i] = max((state1[i-1]*EE, state2[i-1]*ME)) * ET
            state2[i] = max((state1[i-1]*EM, state2[i-1]*MM)) * MT
            i+= 1
    print state1
    print state2
    string = ''
    j = len(seq)-1
    while j >= 0:
        if state1[j] > state2[j]:
               string = string + 'E'
               j -= 1
        else:
               string = string + 'M'
               j -= 1
    print("Predicted output")
    print string
    return max(state1[i-1], state2[i-1])

def viterbi(seq, initE, initM, tMat, cMat):
    return viterbiNum(seq, initE, initM, tMat, cMat) / forwardHMM(seq, initE, initM, tMat, cMat)

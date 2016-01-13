#returns a zero matrix of size x*y
def makeMatrix(x,y):
    mat = [[0]*x]
    county = 0
    while county < y:
        mat.append([0]*x)
        county += 1
    return mat

#prints a matrix by line, minimal format
def printm(matrix):
    count = 0
    while count < len(matrix):
        print matrix[count]
        count += 1
    return

#returns the sum of "boxed area"
def alignSum(matrix,x,y):
    alignArray = []
    countx = x - 1
    while countx >= 0:
        alignArray.append(matrix[countx][y - 1])
        countx -= 1
    county = y - 1
    while county >= 0:
        alignArray.append(matrix[x - 1][county])
        county -= 1
    return sum(alignArray)- matrix[x-1][y-1]

#returns the max of the "boxed area"
def alignMax(matrix,x,y):
    alignArray = []
    countx = x - 1
    while countx >= 0:
        alignArray.append(matrix[countx][y - 1])
        countx -= 1
    county = y - 1
    while county >= 0:
        alignArray.append(matrix[x - 1][county])
        county -= 1
    return max(alignArray)

#matrix initializer for numAligns
def numAlignsMatrix(x,y):
    first = [1]*x
    others = [1]
    others.extend([0]*(x-1))
    county = 1
    mat = [first]
    while county < y:
        mat.append(list(others))
        county += 1
    return mat

#returns the number of alignments possible given two strings
def numAligns(seq1,seq2):
    x = len(seq1)
    y = len(seq2)
    m = numAlignsMatrix(x+1,y+1)
    countx = 1
    county = 1
    while countx <= y:
        while county <= x:
            m[countx][county] = alignSum(m,countx,county)
            county += 1
        countx += 1
        county = 1
    printm(m)
    return m[y][x]

#scoring function for sequence alignment
def score(seq1, seq2, idx1, idx2, match, mismatch):
    if seq1[idx1] == seq2[idx2]:
        return match
    else:
        return mismatch

def evalNeighbours(m,x,y):
    nlist = []
    try:
        nlist.append(m[x-1][y-1])
    except:
        nlist.append(float("-infinity"))
    try:
        nlist.append(m[x-1][y])
    except:
        nlist.append(float("-infinity"))
    try:
        nlist.append(m[x][y-1])
    except:
        nlist.append(float("-infinity"))
    return nlist

def traceback(pwsaMatrix, seq1, seq2):
    m = pwsaMatrix
    string1 = ""
    string2 = ""
    x = len(seq1) - 1
    y = len(seq2) - 1
    while x != 0 and y != 0:
        pos = m[x][y]
        if pos == 0:
            string1 = seq1[x] + string1
            string2 = seq2[y] + string2
            x -= 1
            y -= 1
        if pos < 0:
            string1 = "_"*pos + string1
            count = 0
            while count > pos:
                string2 = seq2[y] + string2
                y -= 1
                count -= 1
        else:
            count = 0
            while count < pos:
                string1 = seq1[x] + string1
                x -= 1
                count += 1
            string2 = "_"*pos + string2
    alignment = [string1, string2]
    return alignment

def tracebackMatrix(x,y):
    first = range(0,x*-1,-1)
    m = [first]
    county = 1
    while county < y:
        other = [county]
        other.extend([0]*(x-1))
        m.append(list(other))
        county += 1
    return m

def PWSA(seq1, seq2, match, mismatch, gap):
    first = range(0,(len(seq1)+1)*(gap),gap)
    mScores =[first]
    #initialize matrix with edges according to penalty
    for idx, val in enumerate(seq2[:]):
        toMult = idx + 1
        other = [toMult*gap]
        other.extend([0]*(len(seq1)))
        mScores.append(list(other))
    #update matrix
    countx = 1
    county = 1
    mScores[1][1] = 1
    while countx <= len(seq2):
        while county <= len(seq1):
            mScores[countx][county] += alignMax(mScores,countx,county) + score(seq1,seq2,county-1,countx-1, match, mismatch)
            county += 1
        countx += 1
        county = 1
##    mTrace = tracebackMatrix(len(seq1),len(seq2))
##    countx = 1
##    county = 1
##    while countx < len(seq2):
##        while county < len(seq1):
##            mTrace[countx][county] = alignSum(mScores,countx,county)
##            county += 1
##        countx += 1
##        county = 1
##    print "Trace Matrix"
##    printm(mTrace)
##    print "------"
    return mScores
#############
m = PWSA("ABCD","ABCD",1,-1,-1)
printm(m)
#print "Diff Length Strings"
#m2 = PWSA("ABCD","ABCDE",1,-1,-2)
#printm(m2)

#function readFile(textfile)
#reads textfile by line and stores them in list 'lines' with each line as a string
#variable 'filename' should be a string of the location of the target text file
def readFile(filename):
    with open(filename) as f:
        lines = f.read().splitlines()
    return lines
#function common_substr(1st string, 2nd string, length of substring)
# taken from http://stackoverflow.com/questions/16447839/common-substring-of-length-k
# finds substrings of length k between two strings and returns them
# only returns one substring per comparison - the first found
def common_substr(a, b, k):
  for substr in (a[i:i+k] for i in range(len(a)-k+1)):
    if substr in b:
      return substr

#function fastKmers(sequences, kmer length, kmers to analyze)
#takes a list of sequences (sequenceList) and runs common_substr on all of them
#   until it finds n kmers and checks for only those kmers in the sequenceList
#found kmers are stored as keys in kdict{}
#every time the kmer is found in a sequence, the value for the key is increased
#the average number of matches per kmer is used as coverage for the estimation of length
def fastKmers(sequenceList, kLength = 17, kmers = 80):
    kdict = {}
    countRow = 0
    countCol = 1
    rowMatched = []
    while countRow < len(sequenceList) -1 and len(kdict.keys()) < kmers:
        while countCol < len(sequenceList) and len(kdict.keys()) < kmers:
            a = common_substr(sequenceList[countRow], sequenceList[countCol], kLength)
            if a:
                kdict[a] = 0
                #print 'Kmer found: ' + a + ' at ' + 'seq' + str(countRow) + ' and seq' + str(countCol)
            countCol += 1
        countRow += 1
        countCol = countRow + 1
    for kmer in kdict:
        for sequence in sequenceList:
            if kmer in sequence:
                kdict[kmer] += 1
    kmerCount = 0.0
    for value in kdict.itervalues():
        kmerCount += value
    #print kmerCount
    coverage = kmerCount / float(len(kdict.keys()))
    #print coverage
    seqLenEstimate = float(len(sequenceList)) * (float(len(sequenceList[0])) - float(kLength) + 1.0) / coverage
    #print kdict
    return seqLenEstimate

####   Main
##
##print 'Estimation of sequence length in \'mySequence.txt\':'
##file1 = readFile('mySequence.txt')
##test1 = fastKmers(file1)
##print test1
##print 'Estimation of sequence length in \'segments.txt\':'
##file2 = readFile('segments.txt')
##test2 = fastKmers(file2)
##print test2
##print 'Estimation of sequence length in \'reads.txt\':'
##file3 = readFile('reads.txt')
##test3 = fastKmers(file3)
##print test3

import re
import random
import operator
from time import process_time


## Functions
def multipleSeedsGibbsSampling(sequences, numSeeds, motifLength, numIter):
    bestScore = 0
    bestMotifs = None
    for i in range(numSeeds):
        results = gibbsSampler(sequences, motifLength, numIter)
        if i == 0:
            bestScore = score(results)
            bestMotifs = results
        else:
            if score(results) < bestScore:
                bestScore = score(results)
                bestMotifs = results
    return bestMotifs


def gibbsSampler(sequences, motifLength, numIter):
    # The number of sequences
    numSeqs = len(sequences)

    # Randomly select motif in each sequence first
    motifs = []
    for seq in sequences:
        i = random.randint(0, len(seq) - motifLength)
        motif = seq[i:i+motifLength]
        motifs.append(motif)
    
    bestMotif = motifs
    bestMotifsScore = score(motifs)

    # Randomly choose a sequence "numIter" times and update motifs or not
    for j in range(numIter):
        i = random.randrange(numSeqs)

        # Compute the profile for the sequences except for "i"
        subsetMotifs = motifs[0:i]+motifs[i+1:numSeqs]
        updatedMotif = profile(subsetMotifs, sequences[i])
        motifs[i] = updatedMotif

        if score(motifs) < bestMotifsScore:
            bestMotifs = motifs
            bestMotifsScore = score(motifs)

    return bestMotifs


def profile(subsetMotifs, sequence):
    # Length of motif
    k = len(subsetMotifs[0])

    # Make position specific score matrix first
    pssm = makePSSM(subsetMotifs)
	
    # Calculate probilities for each k-mer in sequences[i]
    probs = []
    sumProbs = 0.0
    for i in range(len(sequence)-k+1):
        prob = 1.0
        subSeq = sequence[i:i+k]
        for j in range(k):
            if subSeq[j] == 'A':
                prob *= pssm[0][j]
            elif subSeq[j] == 'C':
                prob *= pssm[1][j]
            elif subSeq[j] == 'G':
                prob *= pssm[2][j]
            elif subSeq[j] == 'T':
                prob *= pssm[3][j]
        probs.append(prob)
        sumProbs += prob

    # Randomly select a k-mer
    partialSum = 0.0
    randVal = random.random()
    position = 0
    for i in range(len(sequence)-k+1):
        partialSum += probs[i]
        if partialSum/sumProbs >= randVal:
            position = i
            break

    return sequence[position:position+k]


def makePSSM(subsetMotifs):
    # length of a motif
    k = len(subsetMotifs[0])

    # Initialize pssm by (k x 4)
    pssm = [[0.0 for i in range(k)] for j in range(4)]
    for ithPos in range(k):
        A = C = G = T = 1.0
        for motif in subsetMotifs:
            if motif[ithPos]=='A':
                A += 1.0
            elif motif[ithPos]=='C':
                C += 1.0
            elif motif[ithPos]=='G':
                G += 1.0
            elif motif[ithPos]=='T':
                T += 1.0
	
        total = A + C + G + T
        pssm[0][ithPos] = A/total
        pssm[1][ithPos] = C/total
        pssm[2][ithPos] = G/total
        pssm[3][ithPos] = T/total
    
    return pssm


def score(motifs):
    # length of a motif
    k = len(motifs[0])

    pattern = []
    for ithPos in range(k):
        A = C = G = T = 0
        for motif in motifs:
            if motif[ithPos]=='A':
                A += 1
            elif motif[ithPos]=='C':
                C += 1
            elif motif[ithPos]=='G':
                G += 1
            elif motif[ithPos]=='T':
                T += 1

        if A >= C and A >= G and A >= T:
            pattern.append('A')
        elif C >= G and C >= T:
            pattern.append('C')
        elif G >= T:
            pattern.append('G')
        else:
            pattern.append('T')

    pattern = "".join(pattern)
    score = 0
    for motif in motifs:
        score += hammingDistance(motif, pattern)

    return score


def hammingDistance(str1, str2):
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs


## Main
# Read the sequence text file first
fileQuery = open("/mnt/ssd0/semin/motif-finder/dataset/mm9Gata4MotifCollection.txt", "r")
sequences = []
for line in fileQuery:
    if line[0:3] != '>mm':
        line = line.strip()
        sequences.append(line)
fileQuery.close()

fileAnswer = open("/mnt/ssd0/semin/motif-finder/dataset/mm9Gata4Solutions.txt", "r")
answers = []
for line in fileAnswer:
    if line[0:3] != '>mm':
        line = re.sub('[^A-Z]', '', line)
        answers.append(line)
fileAnswer.close()

# Set the parameters
numSeeds = 20
motifLength = 11
numIter = 1000

# Run gibbsSampler
start_time = process_time()
bestMotifs = multipleSeedsGibbsSampling(sequences, numSeeds, motifLength, numIter)
end_time = process_time()
elapsed_time = end_time - start_time

print("Motif Finding via Gibbs Sampling")
print("System Best Pick      Match/Unmatch      Answer")
for i in range(len(bestMotifs)):
    if bestMotifs[i] == answers[i]:
        print(bestMotifs[i],"          Match             ", answers[i])
    else:
        print(bestMotifs[i],"          Unmatch           ", answer[i])
print("Elapsed Time:", elapsed_time)

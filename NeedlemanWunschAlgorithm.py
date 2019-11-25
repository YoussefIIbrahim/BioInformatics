import sys
import getopt
import argparse
from Bio import SeqIO

GAP_PENALTY = -1
MATCH = 1
MISMATCH = -1
SEQ1 = ''
SEQ2 = ''
MAX_PATHS = 10
MAX_SEQ = 15

seq1List = []
seq2List = []
listOfAlignmentsSeq1 = []
listOfAlignmentsSeq2 = []

finalList = []


def main():

    argv = sys.argv[1:]

    try:
        # Define the getopt parameters
        opts, args = getopt.getopt(argv, 'a:b:c:', ['soperand', 'soperand', 'soperand'])
        # Check if the options' length is 2 (can be enhanced)
        if len(opts) == 0 or len(opts) > 3:
            sys.exit('usage: Needleman-Wunsch.py -a <config.txt path> -b <sequence1.txt path> -c <sequence2.txt path>')

    except getopt.GetoptError:
        sys.exit('usage: Needleman-Wunsch.py -a <config.txt path> -b <sequence1.txt path> -c <sequence2.txt path>')

    ap = argparse.ArgumentParser()

    ap.add_argument("-a", "--apath", required=True, help="<config.txt path>, path to config.txt")

    ap.add_argument("-b", "--bpath", required=True,
                    help="<sequence1.txt path>, path to sequence1.txt")

    ap.add_argument("-c", "--cpath", required=True,
                    help="<sequence2.txt path>, path to sequence2.txt")

    args = vars(ap.parse_args())

    filePath = args['apath']
    sequence1 = args['bpath']
    sequence2 = args['cpath']

    GAP_PENALTY, MATCH, MISMATCH, MAX_PATHS, MAX_SEQ = readFile(str(filePath), 1)

    GAP_PENALTY = int(GAP_PENALTY)
    MATCH = int(MATCH)
    MISMATCH = int(MISMATCH)
    MAX_PATHS = int(MAX_PATHS)
    MAX_SEQ = int(MAX_SEQ)

    # SEQ1 = readFile(str(sequence1), 0)[0]
    # SEQ2 = readFile(str(sequence2), 0)[0]
    SEQ1 = SeqIO.read(sequence1, "fasta")
    SEQ2 = SeqIO.read(sequence2, "fasta")
    try:
        assert len(SEQ1) <= MAX_SEQ and len(SEQ2) <= MAX_SEQ
    except:
        sys.exit("ERROR: MAX_SEQ in config.txt is smaller than actual sequences.")

    output = needleman_wunsch(SEQ1, SEQ2, GAP_PENALTY, MATCH, MISMATCH, MAX_PATHS)
    outputString = 'Score: \n ' + str(output[len(output)-1]) + '\n Alignments: \n'
    for i in range(len(output)):
        if i != len(output)-1:
            for j in range(len(output[i])-1):
                outputString += '(' + output[i][j] + ', ' + output[i][j+1] + ') \n'
    print("File 'output.txt' is saved and the content is:\n", outputString)
    with open("output.txt", 'w') as f:
        f.write(outputString)


def readFile(path, mode):
    configList = []
    try:
        with open(path) as fp:
            line = fp.readline()
            while line:
                val = line.strip()
                if mode:
                    configList.append(val[val.find('= ')+len('= '):])
                else:
                    configList.append(val)
                line = fp.readline()
        return configList
    except:
        sys.exit("ERROR: Couldn't find file. Please make sure you included the correct paths.")


# A function for making a matrix of zeroes
def zeros(rows, cols):
    retval = []
    for x in range(rows):
        retval.append([])
        for y in range(cols):
            retval[-1].append(0)
    return retval


# A function for determining the score between any two bases in alignment
def match_score(alpha, beta, GAP_PENALTY, MATCH, MISMATCH):
    if alpha == beta:
        return MATCH
    elif alpha == '-' or beta == '-':
        return GAP_PENALTY
    else:
        return MISMATCH


def needleman_wunsch(seq1, seq2, GAP_PENALTY, MATCH, MISMATCH, MAX_PATHS):
    n = len(seq1)
    m = len(seq2)

    score = zeros(m + 1, n + 1)

    for i in range(0, m + 1):
        score[i][0] = GAP_PENALTY * i

    for j in range(0, n + 1):
        score[0][j] = GAP_PENALTY * j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + match_score(seq1[j - 1], seq2[i - 1], GAP_PENALTY, MATCH, MISMATCH)
            delete = score[i - 1][j] + GAP_PENALTY
            insert = score[i][j - 1] + GAP_PENALTY

            score[i][j] = max(match, delete, insert)

    i = m
    j = n
    recursionAlignment(seq1, seq2, i, j, score, GAP_PENALTY, MATCH, MISMATCH)

    for i in range(len(listOfAlignmentsSeq1)):
        listOfAlignmentsSeq1[i].reverse()
        listOfAlignmentsSeq2[i].reverse()

    for i in range(len(listOfAlignmentsSeq1)):
        finalList.append((''.join(listOfAlignmentsSeq1[i]), ''.join(listOfAlignmentsSeq2[i])))

    try:
        assert MAX_PATHS >= len(finalList)
    except:
        sys.exit("ERROR: MAX_PATHS in config.txt is smaller than actual value.")
    finalList.append(score[-1][len(score[1])-1])
    return finalList


def recursionAlignment(seq1, seq2, i, j, score, GAP_PENALTY, MATCH, MISMATCH):
    if i == 0 and j == 0:
        listOfAlignmentsSeq1.append(seq1List[:])
        listOfAlignmentsSeq2.append(seq2List[:])
    else:
        score_current = score[i][j]
        score_diagonal = score[i - 1][j - 1]
        score_up = score[i][j - 1]
        score_left = score[i - 1][j]

        if score_current == score_diagonal + match_score(seq1[j - 1], seq2[i - 1], GAP_PENALTY, MATCH, MISMATCH):
            seq1List.append(seq1[j - 1])
            seq2List.append(seq2[i - 1])
            recursionAlignment(seq1, seq2, i-1, j-1, score, GAP_PENALTY, MATCH, MISMATCH)
            seq1List.pop()
            seq2List.pop()

        if score_current == score_up + GAP_PENALTY:
            seq1List.append(seq1[j - 1])
            seq2List.append('-')
            recursionAlignment(seq1, seq2, i, j-1, score, GAP_PENALTY, MATCH, MISMATCH)
            seq1List.pop()
            seq2List.pop()

        if score_current == score_left + GAP_PENALTY:
            seq1List.append('-')
            seq2List.append(seq2[i - 1])
            recursionAlignment(seq1, seq2, i-1, j, score, GAP_PENALTY, MATCH, MISMATCH)
            seq1List.pop()
            seq2List.pop()

main()

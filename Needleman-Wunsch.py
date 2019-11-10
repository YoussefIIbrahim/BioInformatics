gap_penalty = -1
match_award = 1
mismatch_penalty = -1

seq1 = "GCATGCU"
seq2 = "GATTACA"
seq1List = []

seq2List = []

listOfAlignmentsSeq1 = []
listOfAlignmentsSeq2 = []

finalList = []
configList = []


def InitScript():
    try:
        filepath = 'config.txt'
        with open(filepath) as fp:
            line = fp.readline()
            while line:
                val = line.strip()
                configList.append(val[val.find('= ')+len('= '):])
                line = fp.readline()
    except:
        print("ERROR: Couldn't find config.txt file. Please make sure it is in the same directory as the script")

# A function for making a matrix of zeroes
def zeros(rows, cols):
    # Define an empty list
    retval = []
    # Set up the rows of the matrix
    for x in range(rows):
        # For each row, add an empty list
        retval.append([])
        # Set up the columns in each row
        for y in range(cols):
            # Add a zero to each column in each row
            retval[-1].append(0)
    # Return the matrix of zeros
    return retval


# A function for determining the score between any two bases in alignment
def match_score(alpha, beta):
    if alpha == beta:
        return match_award
    elif alpha == '-' or beta == '-':
        return gap_penalty
    else:
        return mismatch_penalty


def needleman_wunsch(seq1, seq2):
    n = len(seq1)
    m = len(seq2)

    score = zeros(m + 1, n + 1)

    for i in range(0, m + 1):
        score[i][0] = gap_penalty * i

    for j in range(0, n + 1):
        score[0][j] = gap_penalty * j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + match_score(seq1[j - 1], seq2[i - 1])
            delete = score[i - 1][j] + gap_penalty
            insert = score[i][j - 1] + gap_penalty

            score[i][j] = max(match, delete, insert)

    i = m
    j = n

    recursionAlignment(seq1, seq2, i, j, score)

    for i in range(len(listOfAlignmentsSeq1)):
        listOfAlignmentsSeq1[i].reverse()
        listOfAlignmentsSeq2[i].reverse()

    for i in range(len(listOfAlignmentsSeq1)):
        finalList.append((''.join(listOfAlignmentsSeq1[i]), ''.join(listOfAlignmentsSeq2[i])))

    print(finalList)
    return finalList


def recursionAlignment(seq1, seq2, i, j, score):

    if not(i > 0 and j > 0):
        listOfAlignmentsSeq1.append(seq1List[:])
        listOfAlignmentsSeq2.append(seq2List[:])
    else:
        score_current = score[i][j]
        score_diagonal = score[i - 1][j - 1]
        score_up = score[i][j - 1]
        score_left = score[i - 1][j]

        if score_current == score_diagonal + match_score(seq1[j - 1], seq2[i - 1]):
            seq1List.append(seq1[j - 1])
            seq2List.append(seq2[i - 1])
            recursionAlignment(seq1, seq2, i-1, j-1, score)
            seq1List.pop()
            seq2List.pop()

        if score_current == score_up + gap_penalty:
            seq1List.append(seq1[j - 1])
            seq2List.append('-')
            recursionAlignment(seq1, seq2, i, j-1, score)
            seq1List.pop()
            seq2List.pop()

        if score_current == score_left + gap_penalty:
            seq1List.append('-')
            seq2List.append(seq2[i - 1])
            recursionAlignment(seq1, seq2, i-1, j, score)
            seq1List.pop()
            seq2List.pop()

InitScript()
needleman_wunsch(seq1, seq2)

# print(output1 + "\n" + output2)

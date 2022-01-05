import blosum62
blosum62 = blosum62.blosum62

#Local DNA
def localDNA(seq_1, seq_2):
    seq_1_len = len(seq_1)
    seq_2_len = len(seq_2)
    negativeScore = 0
    matchScore = 1
    misMatchScore = -2
    gapScore = -1
    maxScore = 0
    maxScorePosition = 0
    alignmentMatrix = [[0 for seq_1 in range(seq_1_len + 1)] for seq_2 in range(seq_2_len + 1)]

    for row in range(1, seq_2_len + 1, 1):
        for col in range(1, seq_1_len + 1, 1):
            if seq_1[col - 1] == seq_2[row - 1]:
                alignmentMatrix[row][col] = alignmentMatrix[row - 1][col - 1] + matchScore
            else:
                alignmentMatrix[row][col] = max(alignmentMatrix[row - 1][col] + gapScore,
                                                alignmentMatrix[row][col - 1] + gapScore,
                                                alignmentMatrix[row - 1][col - 1] + misMatchScore,
                                                negativeScore)
            if alignmentMatrix[row][col] > maxScore:
                maxScore = alignmentMatrix[row][col]
                maxScorePosition = (row, col)
    print(alignmentMatrix)
    # tracing back
    Al_seq_1 = ""
    Al_seq_2 = ""
    match = ""
    row = maxScorePosition[0]
    col = maxScorePosition[1]
    while alignmentMatrix[row][col] != 0:

        nextCell = max(alignmentMatrix[row - 1][col],
                       alignmentMatrix[row][col - 1],
                       alignmentMatrix[row - 1][col - 1])
        # go to left cell
        if alignmentMatrix[row - 1][col - 1] == nextCell:
            Al_seq_2 += seq_2[row - 1]
            Al_seq_1 += seq_1[col - 1]
            match += "|"
            row -= 1
            col -= 1
        elif alignmentMatrix[row][col - 1] == nextCell:
            Al_seq_1 += seq_1[col - 1]
            Al_seq_2 += "_"
            match += " "
            col -= 1
        elif alignmentMatrix[row - 1][col] == nextCell:
            Al_seq_1 += "_"
            Al_seq_2 += seq_2[row - 1]
            match += " "
            row -= 1

    Al_seq_1 = Al_seq_1[::-1]
    Al_seq_2 = Al_seq_2[::-1]
    match = match[::-1]
    print(Al_seq_1)
    print(match)
    print(Al_seq_2)



###############################################################

#Local Protein

# get protein score
def getScore(firstProtein, secondProtein):
    if(firstProtein,secondProtein) in blosum62:
        score = blosum62[(firstProtein, secondProtein)]
    else:
        score = blosum62[(secondProtein, firstProtein)]
    
    return score

# local protein alignment
def getAlignedSeqs(seq_1, seq_2):
    seq_1_len = len(seq_1)
    seq_2_len = len(seq_2)
    maxScore = 0
    maxScorePosition = 0

    # Initialize alignment matrix
    alignmentMatrix = [[0 for seq_1 in range(seq_1_len + 1)] for seq_2 in range(seq_2_len + 1)]

    # Fill alignment matrix with protein score
    for row in range(1, seq_2_len + 1, 1):
        for col in range(1, seq_1_len + 1, 1):
            score = getScore(seq_1[col - 1], seq_2[row - 1])
            if seq_1[col - 1] == seq_2[row - 1]:
                alignmentMatrix[row][col] = alignmentMatrix[row - 1][col - 1] + score
            else:
                alignmentMatrix[row][col] = max(alignmentMatrix[row - 1][col] + score,
                                                alignmentMatrix[row][col - 1] + score,
                                                alignmentMatrix[row - 1][col - 1] + score, 0)
            if alignmentMatrix[row][col] > maxScore:
                maxScore = alignmentMatrix[row][col]
                maxScorePosition = (row, col)

    #print(alignmentMatrix)

    # tracing back
    Al_seq_1 = ""
    Al_seq_2 = ""
    match = ""
    row = maxScorePosition[0]
    col = maxScorePosition[1]
    while alignmentMatrix[row][col] != 0:

        nextCell = max(alignmentMatrix[row - 1][col],
                       alignmentMatrix[row][col - 1],
                       alignmentMatrix[row - 1][col - 1])
        if alignmentMatrix[row - 1][col - 1] == nextCell:
            Al_seq_2 += seq_2[row - 1]
            Al_seq_1 += seq_1[col - 1]
            match += "|"
            row -= 1
            col -= 1
        elif alignmentMatrix[row][col - 1] == nextCell:
            Al_seq_1 += seq_1[col - 1]
            Al_seq_2 += "_"
            match += " "
            col -= 1
        elif alignmentMatrix[row - 1][col] == nextCell:
            Al_seq_1 += "_"
            Al_seq_2 += seq_2[row - 1]
            match += " "
            row -= 1

    Al_seq_1 = Al_seq_1[::-1]
    Al_seq_2 = Al_seq_2[::-1]
    match = match[::-1]
    print(Al_seq_1)
    print(match)
    print(Al_seq_2)



#Start Program
while(True):

    print("1: Local Allignment DNA")
    print("2: Local Allignment Protein")
    print("3: Close")

    number = int(input("Choose : "))

    if number == 1:
        seq1 = input("Enter the first DNA sequence: ")
        seq2 = input("Enter the second DNA sequence:")
        localDNA(seq1, seq2)
        print("------------------------------------------------")
    elif number == 2:
        seq1 = input("Enter the first protein sequence: ")
        seq2 = input("Enter the second protein sequence:")
        getAlignedSeqs(seq1, seq2)
        print("------------------------------------------------")
    elif number == 3:
        break
    else:
        print("Choose valid Number")
        print("------------------------------------------------")

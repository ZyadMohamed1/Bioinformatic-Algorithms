import blosum62
from collections import defaultdict

b62 = blosum62.blosum62
HSP_withScore = defaultdict(dict)
SubstringsDict = {}
SeqID = {}


# step 1
def cleanSeq(querySeq):
    pairsList = []
    cleaned = ""
    for i in range(0, len(querySeq), 2):
        pairsList.append(querySeq[i: i + 2])

    for i in range(0, len(pairsList) - 1):
        if pairsList[i] == pairsList[i + 1]:
            pairsList[i] = "XX"
            pairsList[i + 1] = "XX"

    for z in range(0, len(pairsList)):
        cleaned = cleaned + pairsList[z]

    return cleaned


# step 2
def getWords(querySeq, wordLength):
    seqLength = len(querySeq)
    index = 0
    wordList = []
    while index < seqLength and index + wordLength <= seqLength:
        wordList.append(querySeq[index: index + wordLength])
        index += 1

    return wordList


# step 3 and 4
def getSeeds(wordList, threshHold):

    letters = blosum62.letters
    seeds = []

    # loop for the wordsList
    for i in range(len(wordList)):
        word = wordList[i]
        # Loop for the letters of the word
        for j in range(len(word)):
            letter = 0
            # Loop for all Letters to be changed
            for k in range(len(letters)):
                score = 0
                wordChange = list(word)
                wordChange[j] = letters[letter]
                new_word = "".join(wordChange)

                # Summing the blosum62
                for l in range(len(word)):
                    if (new_word[l], word[l]) in b62:
                        score = b62[(new_word[l], word[l])] + score
                    else:
                        score = b62[(word[l], new_word[l])] + score

                if score >= threshHold:
                    seeds.append(new_word)
                    SubstringsDict[new_word] = word

                letter = letter + 1

    return seeds


# Search for hit words
def searchInFile(seeds):
    hit_words = []
    for i in range(len(seeds)):
        with open('database.txt') as file:
            if seeds[i] in file.read():
                hit_words.append(seeds[i])

    return hit_words


# get list of database sequences
def getDataBase():
    file = open("database.txt", "r")
    lines = file.readlines()
    seqNum = 1
    DBSequences = []
    seq = ""
    for i in lines:
        if i.__contains__(">") or i == lines[len(lines) - 1]:
            if seq != "":
                DBSequences.append(seq)
            seq = ""
            seqNum += 1
        else:
            seq = seq + i[0:len(i) - 1]
    return DBSequences


# Word extension
def wordExtension(hitWords, dataBaseList, querySequence, HSP):
    ID = 0
    for sequences in range(len(dataBaseList)):
        SeqID[ID] = 0
        dataBaseSequence = dataBaseList[sequences]

        for i in range(len(hitWords)):
            rep = 0
            if hitWords[i] in dataBaseSequence:
                
                while True:

                    score = 0
                    bestScore = 0
                    startIndexDatabase = dataBaseSequence.find(hitWords[i], rep)
                    lastIndexDatabase = startIndexDatabase + (len(hitWords[i]) - 1)
                    startIndexQuery = querySequence.find(SubstringsDict.get(hitWords[i]))
                    lastIndexQuery = startIndexQuery + len(SubstringsDict.get(hitWords[i])) - 1
                    rep = startIndexDatabase
                    word = hitWords[i]
                    if rep == -1:
                        break
                    else:
                        rep = rep + 1
                        for j in range(len(word)):
                            score = score + b62[(word[j], word[j])]

                    if score >= HSP:
                        bestScore = score
                        # loop to extend left and right
                        while (startIndexQuery > 0 and startIndexDatabase > 0) \
                                or (lastIndexQuery < len(querySequence) - 1
                                    and lastIndexDatabase < len(dataBaseSequence) - 1):

                            tmp = score
                            tmpSequence = word

                            # Extend left
                            if startIndexQuery > 0 and startIndexDatabase > 0:
                                if (querySequence[startIndexQuery - 1], dataBaseSequence[startIndexDatabase - 1]) in b62:
                                    score = b62[(querySequence[startIndexQuery - 1],
                                                dataBaseSequence[startIndexDatabase - 1])] + score
                                else:
                                    score = b62[(dataBaseSequence[startIndexDatabase - 1],
                                                querySequence[startIndexQuery - 1])] + score

                                word = querySequence[startIndexQuery - 1] + word
                                startIndexQuery = startIndexQuery - 1
                                startIndexDatabase = startIndexDatabase - 1

                            # Extend right
                            if lastIndexQuery < len(querySequence) - 1 and lastIndexDatabase < len(dataBaseSequence) - 1:
                                if (querySequence[lastIndexQuery + 1], dataBaseSequence[lastIndexDatabase + 1]) in b62:
                                    score = b62[(querySequence[lastIndexQuery + 1],
                                                dataBaseSequence[lastIndexDatabase + 1])] + score
                                else:
                                    score = b62[(dataBaseSequence[lastIndexDatabase + 1],
                                                querySequence[lastIndexQuery + 1])] + score

                                word = word + querySequence[lastIndexQuery + 1]
                                lastIndexQuery = lastIndexQuery + 1
                                lastIndexDatabase = lastIndexDatabase + 1

                            # check score
                            if score >= HSP:
                                if score >= bestScore:
                                    bestScore = score

                                if (bestScore - score) >= 3:
                                    bestScore = tmp
                                    word = tmpSequence
                                    break
                            else:
                                bestScore = score
                                break
                        #
                        if ID in HSP_withScore:
                            if word in HSP_withScore[ID]:
                                HSP_withScore[ID][word] = HSP_withScore[ID][word] + bestScore
                            else:
                                HSP_withScore[ID][word] = bestScore
                                # Add Non repeated sum 1
                                # #SeqID[ID] = SeqID[ID] + bestScore
                        else:
                            HSP_withScore[ID][word] = bestScore
                            # Add Non repeated sum 2
                            # SeqID[ID] = SeqID[ID] + bestScore
                        # If i want to add repeated sum
                        SeqID[ID] = SeqID[ID] + bestScore
        ID = ID + 1


# Input
querySequence = input("Enter Query Sequence :")
words = int(input("Enter substring length : "))
ThreshHold = int(input("Enter T : "))
HSP = int(input("Enter HSP ThreshHold: "))

# step #1 Removing Low complexity regions
cleanedSeq = cleanSeq(querySequence)
# Step #2 get words list from the query sequence
wordsList = getWords(cleanedSeq, words)
# Steps #3 and #4 find the neighborhood words and get the seeds
Seed = getSeeds(wordsList, ThreshHold)
Seed = list(dict.fromkeys(Seed))
# Step #5 scanning the database
HitWords = searchInFile(Seed)
# Step #6 and #7
Database = getDataBase()
wordExtension(HitWords, Database, querySequence, HSP)

# Output
print(HitWords)
print("Each HSE: ")
print(HSP_withScore)
print("Sequences Total HSE: ")
print(SeqID)
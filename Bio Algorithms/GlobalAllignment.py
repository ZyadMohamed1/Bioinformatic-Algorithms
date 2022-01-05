import blosum62
blosum62 = blosum62.blosum62

#Global DNA and Protein
def AllignDis(Fseq, Sseq, match , mismatch , gap_penalty ):
    first=len(Fseq)
    second=len(Sseq)

    matrix = [[0 for x in range(first+1)] for y in range(second+1)] 
    for i in range (1,second+1):
        matrix[i][0]=matrix[i-1][0] + gap_penalty
    for j in range (1,first+1):
        matrix[0][j] =matrix[0][j-1] + gap_penalty

    if Type == "D":
        for i in range(1,second+1):
            for j in range(1,first+1):
                if Fseq[j-1] == Sseq[i-1]:
                 score = match
                else :
                    score = mismatch
                matrix[i][j] = max(matrix[i-1][j]+gap_penalty, matrix[i][j-1]+gap_penalty,matrix[i-1][j-1]+score)
    elif Type == "P":
        for i in range(1,second+1):
            for j in range(1,first+1):
                if(Sseq[i-1], Fseq[j-1]) in blosum62:
                    score= blosum62[(Sseq[i-1], Fseq[j-1])]
                else:
                    score= blosum62[(Fseq[j-1], Sseq[i-1])]
                matrix[i][j] =  max(matrix[i-1][j]+gap_penalty, matrix[i][j-1]+gap_penalty,matrix[i-1][j-1]+score)   
    print(matrix[second][first])


    FAllign =""
    SecAllign=""
    MatchAllign=""

    first=len(Fseq)
    second=len(Sseq)

    i = second
    j = first
    while i > 0 and j > 0:
        up = matrix[i-1][j]+gap_penalty
        left=matrix[i][j-1]+gap_penalty
        if Type == "D":
            if Fseq[j-1]==Sseq[i-1]:
                score = match
            else :
                score = mismatch
        elif Type == "P":
            if(Sseq[i-1], Fseq[j-1]) in blosum62:
                score= blosum62[(Sseq[i-1], Fseq[j-1])]
            else:
                score= blosum62[(Fseq[j-1], Sseq[i-1])]

                
        diagonal = matrix[i-1][j-1]+score
        if matrix[i][j]==diagonal:
            FAllign+=Fseq[j-1]
            SecAllign+=Sseq[i-1]
            if score == match:
                MatchAllign+="|"
            else:
                MatchAllign+=" "
            i-=1
            j-=1
        elif matrix[i][j] == up:
            FAllign+="-"
            MatchAllign+=" "
            SecAllign+=Sseq[i-1]
            i-=1
        else:
            FAllign+=Fseq[j-1]
            SecAllign+="-"
            MatchAllign+=" "
            j-=1

    while (i>0):
        FAllign+="-"
        SecAllign+=Sseq[i-1]
        MatchAllign+=" "
        i-=1
    while(j>0):
        FAllign+=Fseq[j-1]
        SecAllign+="-"
        MatchAllign+=" "
        j-=1
        
    FAllign=FAllign[::-1]
    MatchAllign=MatchAllign[::-1]
    SecAllign=SecAllign[::-1]
    
    print (FAllign)
    print(MatchAllign)
    print(SecAllign)

while(True):

    number = int(input("Choose : "))

    if number == 1:
        Type = "D"
        Seq1 = input("Enter First Sequence : ")
        Seq2 = input("Enter Second Sequence : ")
        match = 1
        mismatch = -2
        gap = -1
        AllignDis(Seq1, Seq2, match , mismatch , gap)
        print("------------------------------------------------")
    elif number == 2:
        Type = "P"
        Seq1 = input("Enter First Sequence : ")
        Seq2 = input("Enter Second Sequence : ")
        gap = int(input("Enter The gap Value: "))
        AllignDis(Seq1, Seq2, 0 , 0 , gap)
        print("------------------------------------------------")

    elif number == 3:
        break
    else:
        print("Choose valid Number")
        print("------------------------------------------------")


    
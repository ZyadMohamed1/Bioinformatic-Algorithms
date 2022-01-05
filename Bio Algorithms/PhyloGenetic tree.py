# this fuction gets matrix from file input

def get_matrix():
    File = open('database.txt', 'r')
    matrix = []
    matrix = [line.split() for line in File]

    print(matrix)
    return matrix


# makes label for the matrix to print clustring form

def make_label(Len_matrix):
    seq_list = []
    x = ord("A")
    seq_list.append("A")
    for i in range(1, Len_matrix):
        c = x + i
        seq_list.append(chr(c))
    return seq_list


# splits matrix by the diagonal the form upgma on it

def split_matrix(matrix):
    n = len(matrix)
    split_matrix = [[0 for i in range(0, n)] for j in range(0, n)]
    for i in range(n):
        for j in range(n):
            if (j <= i):
                split_matrix[i][j] = matrix[i][j]
    # print(f)
    for i in range(0, len(split_matrix)):
        for j in range(0, len(split_matrix)):
            split_matrix[i][j] = int(split_matrix[i][j])
    return split_matrix


# search for the minimum value and get its row index and column index

def get_min(matrix, seqList):
    min = 10000000
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if matrix[i][j] < min and matrix[i][j] != 0:
                minI = i
                minJ = j
                min = matrix[i][j]
    print(minI)
    print(minJ)
    merge_matrix(matrix, minI, minJ, seqList)


# preform upgma idea, calculation. Preform and print cluster form

def merge_matrix(matrix, minI, minJ, seq_labels):
    if minJ < minI:
        minI, minJ = minJ, minI
    m = minI
    n = minJ
    for i in range(0, len(matrix)):  # 0 to 5 1

        if (i < minI):

            matrix[minI][i] = (matrix[minI][i] + matrix[minJ][i]) / 2

        elif i > minI and i <= minJ:

            matrix[i][minI] = (matrix[i][minI] + matrix[minJ][i]) / 2

        elif i > minJ and i <= len(matrix):

            matrix[i][minI] = (matrix[i][minI] + matrix[i][minJ]) / 2
            del matrix[i][minJ]

    del matrix[minJ]
    print(matrix)

    seq_labels[m] = "(" + seq_labels[m] + "," + seq_labels[n] + ")"
    # b.append(seqList[m])
    print(seq_labels[m])
    del seq_labels[n]
    return seq_labels

#Start Program
def main():
    matrix = get_matrix()
    len_matrix = len(matrix)
    seqList = make_label(len_matrix)
    matrix = split_matrix(matrix)
    print(matrix)

    while len(matrix) > 1:
        get_min(matrix, seqList)


main()
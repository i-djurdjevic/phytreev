import copy

dStarMatrices = []
dMatrices = []
limbLengthTouples = []
minElements = []
distanceLabels = []
stepDistanceCalculations = []

def parse_steps(steps):
    stepsArray = []
    i = 0
    for elem in steps:
        tmp = (minElements[i],) + elem
        i += 1
        stepsArray.append(tmp)
    return stepsArray

def IfShouldAddCalculation(map, key1, key2):
    if len(map.keys()) != 0 and (key1 not in map) and (key2 not in map):
        return True
    elif len(map.keys()) == 0:
        return True
    else:
        return False

def neighborJoining(matrix, n, labels):
    dMatrices.append(matrix)
    if n == 2:
        edge = {}
        weight = {}
        for i in matrix:
            for j in matrix[i]:
                if j == i:
                    continue
                edge[i] = [j]
                edge[j] = [i]
                weight[(i, j)] = matrix[i][j]
                weight[(j, i)] = matrix[i][j]
        return (edge, weight)
    # Construct neighbor-joining matrix D* from D.
    totalDistance = {}
    for i in matrix:
        totalDistance[i] = 0
        for j in matrix[i]:
            totalDistance[i] += matrix[i][j]
    MatrixD = {}
    dStarCalculations = []
    existingCalculations = {}
    for i in matrix:
        MatrixD[i] = {}
        for j in matrix[i]:
            if i == j:
                MatrixD[i][j] = 0
            else:
                totalDistance[i] = round(totalDistance[i], 2)
                totalDistance[j] = round(totalDistance[j], 2)
                MatrixD[i][j] = round((n - 2) * matrix[i][j] - totalDistance[i] - totalDistance[j], 2)
                # distanceLabels.append((labels[i], labels[j], MatrixD[i][j]))
                dStarCalculation = "(n - 2) * D(" + labels[i] + "," + labels[j] + ")" + " - totalDistance(" + labels[i] + ") - totalDistance(" + labels[j] + ") = " + str(n - 2) + " * " + str(matrix[i][j]) + " - " + str(totalDistance[i]) + " - " + str(totalDistance[j]) + " = " + str(MatrixD[i][j])
                if IfShouldAddCalculation(existingCalculations, "D(" + labels[i] + "," + labels[j] + ")", "D(" + labels[j] + "," + labels[i] + ")"):
                    dStarCalculations.append(dStarCalculation)
                existingCalculations["D(" + labels[i] + "," + labels[j] + ")"] = True
                existingCalculations["D(" + labels[j] + "," + labels[i] + ")"] = True


    if dStarCalculations != []:
        stepDistanceCalculations.append(dStarCalculations)

    step = None
    step = copy.deepcopy(MatrixD)
    dStarMatrices.append(step)
    # Find a minimum element D*i,j of D*.
    mini = float("inf")
    _i = _j = None
    for i in MatrixD:
        for j in MatrixD[i]:
            if MatrixD[i][j] < mini:
                mini = MatrixD[i][j]
                _i, _j = i, j
    minElements.append(mini)
    # Compute delta(i,j) = (TotalDistanceD(i) –TotalDistanceD(j)) / (n – 2).
    delta = (totalDistance[_i] - totalDistance[_j]) / (n - 2.0)
    # Set LimbLength(i) equal to 1/2(delta_i,j + delta_i,j) and LimbLength(j) equal to 1/2(Di,j – delta_i,j).
    limbLength_i = (matrix[_i][_j] + delta) / 2.0
    limbLength_j = (matrix[_i][_j] - delta) / 2.0
    limbLengthTouples.append((limbLength_i, limbLength_j))
    # Form a matrix D' by removing i-th and j-th row/column from D and adding an m-th row/column
    # such that for any k, Dk,m = (Di,k + Dj,k – Di,j) / 2.
    rows = list(matrix)
    m = max(matrix) + 1
    matrix[m] = {}
    matrix[m][m] = 0.0
    for k in rows:
        if k != _i and k != _j:
            matrix[m][k] = (matrix[_i][k] + matrix[_j][k] - matrix[_i][_j]) / 2.0
            matrix[k][m] = matrix[m][k]
            del matrix[k][_i]
            del matrix[k][_j]

    innerLabel = "(" + labels[_i] + "," + labels[_j] + ")"
    labels.append(innerLabel)

    del matrix[_i]
    del matrix[_j]
    # Apply NeighborJoining to D' to obtain Tree(D').
    edge, weight = neighborJoining(matrix, n - 1, labels)

    if len(distanceLabels) == 0:
        distanceLabels.append(labels)
    # Reattach limbs of i and j to obtain Tree(D).
    edge[_i] = [m]
    edge[m].append(_i)
    edge[_j] = [m]
    edge[m].append(_j)
    weight[(_i, m)] = limbLength_i
    weight[(m, _i)] = limbLength_i
    weight[(_j, m)] = limbLength_j
    weight[(m, _j)] = limbLength_j
    return edge, weight

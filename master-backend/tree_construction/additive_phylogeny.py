import copy

import numpy as np

limbLengths = []
limbLengthCaluclations = []
nodes = []
dBalds = []
dTrims = []
minLimbs = []

def limbLength(matrix, j):
    """
    Compute the Limb Length by
    limb = min_ik(d_ij + d_jk - d_ik)/2
    """
    stepLimbLengthCalculations = []
    limb = np.inf
    mask = [i for i in range(matrix.shape[0]) if i != j]
    for k in range(matrix.shape[0]):
        if k == j:
            continue
        limb = min(limb, min(matrix[mask, j] - matrix[mask, k] + matrix[j, k]))
        minLimb = matrix[j, k]
        if len(stepLimbLengthCalculations) == 0:
            length = matrix[mask, j].size
            for i in range(length):
                limbLengthCalculationString = "(" + str(matrix[mask, j][i]) + "-" + str(matrix[mask, k][i]) + "+" + str(
                    matrix[j, k]) + ") / 2 = " + str((matrix[mask, j][i] - matrix[mask, k][i] + matrix[j, k]) / 2)
                stepLimbLengthCalculations.append(limbLengthCalculationString)
            limbLengthCaluclations.append(stepLimbLengthCalculations)
        minLimbs.append(minLimb)
    limbLengths.append(limb // 2)
    nodes.append(j)
    return limb // 2


def find(matrix):
    """
    find i, k s.t. Di,k = Di,n + Dn,k
    """
    for k in range(matrix.shape[0] - 1):
        arr = matrix[k] - matrix[-1]
        index = np.where(arr == matrix[k, -1])
        if len(index[0]) > 0:
            return (index[0][0], k)
    return None


def nearest(edge, weight, x, i, k):
    """
    find the nearest two nodes on path i -> k
    to insert new node, BFS
    """
    queue = [[i]]
    visited = set([i])
    findPath = []
    while len(queue) > 0:
        path = queue.pop()
        node = path[-1]
        visited.add(node)
        if node == k:
            findPath = path
            break
        for next_node in edge[node]:
            if next_node not in visited:
                queue.append(path + [next_node])

    # distance
    dist = 0
    for k in range(len(findPath) - 1):
        i, j = findPath[k], findPath[k + 1]
        if dist + weight[(i, j)] > x:
            return (i, j, x - dist, dist + weight[(i, j)] - x)
        dist += weight[(i, j)]


def additivePhylogeny(matrix, n, inner_n, labels):
    """
    finds the simple tree fitting an n x n
    additive distance matrix D.

    Tree: Edge(u, v) as [... u-th[v, ...]..]
    Weight(u, v) as {(u, v): w}
    inner_n, count the index of inner node
    """
    if n == 2:
        edge = {}
        edge[0] = [1]
        edge[1] = [0]
        weight = {}
        weight[(0, 1)] = matrix[0, 1]
        weight[(1, 0)] = matrix[0, 1]
        return (edge, weight, inner_n)

    limb = limbLength(matrix, n - 1)
    matrix[:-1, -1] -= limb
    matrix[-1, :-1] -= limb

    baldStepMatrix = None
    baldStepMatrix = copy.deepcopy(matrix)
    trimmedStepMatrix = None
    trimmedStepMatrix = copy.deepcopy(baldStepMatrix[:-1, :-1])

    lenDiff = len(baldStepMatrix) - len(labels)
    if lenDiff != 0:
        labels.remove(labels[lenDiff])
    result = np.vstack((labels, baldStepMatrix))
    labels.insert(0, ' ')
    result = np.column_stack((labels, result))
    dBalds.append(result)
    labels.remove(' ')

    labels.remove(labels[lenDiff - 1])
    result = np.vstack((labels, trimmedStepMatrix))
    labels.insert(0, ' ')
    result = np.column_stack((labels, result))
    dTrims.append(result)
    labels.remove(' ')

    # i, n, k three leaves such that Di,k = Di,n + Dn,k
    i, k = find(matrix)
    x = matrix[i, -1]

    # remove row n and column n from D
    edge, weight, inner_n = additivePhylogeny(matrix[:-1, :-1], n - 1, inner_n, labels)
    # the (potentially new) node in T at distance x from i on the path between i and k

    # find the insert node
    i_near, k_near, i_x, n_x = nearest(edge, weight, x, i, k)
    limbLengths.append(i_x)
    # step.limbLength = i_x
    new_node = i_near

    # check if we need to create a new node
    if i_x != 0:
        new_node = inner_n
        inner_n += 1
        # insert between i, k
        edge[i_near].remove(k_near)
        edge[k_near].remove(i_near)
        edge[i_near].append(new_node)
        edge[k_near].append(new_node)
        edge[new_node] = [i_near, k_near]

        weight[(new_node, i_near)] = i_x
        weight[(i_near, new_node)] = i_x
        weight[(new_node, k_near)] = n_x
        weight[(k_near, new_node)] = n_x
        del weight[(i_near, k_near)]
        del weight[(k_near, i_near)]

    # add leaf n back to T by creating a limb (v, n) of length limbLength
    edge[new_node].append(n - 1)
    edge[n - 1] = [new_node]
    weight[(n - 1, new_node)] = limb
    weight[(new_node, n - 1)] = limb
    return (edge, weight, inner_n)

from collections import deque
def depth_first_search(adjacency_matrix, source_node):
    open = deque([source_node])
    parents = {source_node: None}
    result = ""

    while open:
        current_node = open.pop()

        for child_node in adjacency_matrix[current_node]:
            if child_node not in parents.keys():
                parents[child_node] = current_node
                result += str(current_node) + "->" + str(child_node) + "\n"
                open.append(child_node)

    return result

def printOut(edge, weight):
    for i in sorted(edge):
        for j in sorted(edge[i]):
            print(str(i) + "->" + str(j) + ":" + str(weight[(i, j)]))
    # depth_first_search(edge, 0)

class Node():
    def __init__(self, data):
        self.data = data
        self.left = None
        self.right = None
        self.parent = None
        self.vals = [float('inf'), float('inf'), float('inf'), float('inf')]
        self.string = ''
        self.min_par = 0

def toGraph(tree_data):
    graph = {}
    # for input:
    data = tree_data.split('\n')
    leaf = data[0].split("->")[1]
    # get length of alignment for file:
    # leaf = data[0].split("->")[1][:-1]
    li = 0
    for i in range(len(data)):
        # for file:
        # value = data[i][:-1]
        # for input:
        if data[i] != '':
            value = data[i]
            vals = value.split('->')
            vals[0] = int(vals[0])
            if vals[1].isdigit():
                vals[1] = int(vals[1])
            if vals[0] not in graph:
                li = max(li, vals[0])
                graph[vals[0]] = [vals[1]]
            else:
                graph[vals[0]].append(vals[1])

    keys = list(graph.keys())
    nodes = {}
    # nodes[li] = Node(li)
    for n in keys:
        if n not in nodes.keys():
            nodes[n] = Node(n)
            nodes[n].string = str(n)
        for val in graph[n]:
            node = Node(val)
            # if not str(val).isdigit():
            node.string = str(val)
            node.parent = nodes[n]
            if nodes[n].left != None:
                nodes[n].right = node
            else:
                nodes[n].left = node
            nodes[val] = node

    # data = [li]
    # # data_keys = [li]
    # while data != []:
    #     v = data.pop(0)
    #     if v in graph:
    #         data += graph[v]
    #         val = graph[v]
    #         left_node = nodes[val[0]]
    #         left_node.parent = nodes[v]
    #         nodes[v].left = left_node
    #         right_node = nodes[val[1]]
    #         right_node.parent = nodes[v]
    #         nodes[v].right = right_node

    return nodes[keys[0]]

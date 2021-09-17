import utils.newick_utils as nf
class Node():
    def __init__(self, data):
        self.data = data
        self.left = None
        self.right = None
        self.parent = None
        self.vals = [float('inf'), float('inf'), float('inf'), float('inf')]
        self.string = ''
        self.min_par = 0
        self.valsCalculations = ["", "", "", ""]

def post_order(root, index):
    if root:
        # First recur on left child
        post_order(root.left, index)
        # the recur on right child
        post_order(root.right, index)
        if root.left == None:
            print(index)
            print(root.string)
            if root.string[index] == 'A':
                root.vals[0] = 0
            elif root.string[index] == 'C':
                root.vals[1] = 0
            elif root.string[index] == 'G':
                root.vals[2] = 0
            elif root.string[index] == 'T':
                root.vals[3] = 0
        else:
            for i in range(4):
                leftvals = [0, 0, 0, 0]
                rightvals = [0, 0, 0, 0]
                for j in range(4):
                    leftvals[j] = root.left.vals[j] + cost[i][j]
                    rightvals[j] = root.right.vals[j] + cost[i][j]
                root.vals[i] = min(leftvals) + min(rightvals)
                root.valsCalculations[i] += str(min(leftvals)) + " + " + str(min(rightvals)) + " = " + str(root.vals[i])

def level_order(root):
        current = [root]
        root.string += chars[root.vals.index(min(root.vals))]
        root.min_par += min(root.vals)
        while current:
            next_level = []
            for node in current:
                if node.left.left != None:
                    ind = rev_chars[node.string[-1]]
                    v = node.vals[ind]
                    if node.left and node.right:
                        for i in range(4):
                            set = False
                            for j in range(4):
                                if node.left.vals[i] + node.right.vals[j] + cost[ind][i] + cost[ind][j] == v:
                                    node.left.string += chars[i]
                                    node.right.string += chars[j]
                                    set = True
                                    break
                            if set == True:
                                break
                if node.left.left:
                    next_level.append(node.left)
                if node.right.left:
                    next_level.append(node.right)
            current = next_level

def hamming_distance(s1, s2):
    count = 0
    for c1, c2 in zip(s1, s2):
        if c1 != c2:
            count += 1
    return str(count)

def small_parsimony(root, index):
    post_order(root, index)
    level_order(root)

def small_parsimony_wrapper(tree_data):
    graph = {}
    # for input:
    data = tree_data.split('\n')
    # za jedan karakter:
    # leaf = data[0].split("->")[1][0]
    leaf = data[0].split("->")[1]
    # get length of alignment for file:
    # leaf = data[0].split("->")[1][:-1]
    li = 0
    for i in range(len(data)):
        # for file:
        # value = data[i][:-1]
        # for input:
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
    nodes[li] = Node(li)
    root = nodes[li]
    for n in keys:
        for val in graph[n]:
            node = Node(val)
            if not str(val).isdigit():
                node.string = val
            nodes[val] = node

    data = [li]
    while data:
        v = data.pop(0)
        if v in graph:
            data += graph[v]
            val = graph[v]
            left_node = nodes[val[0]]
            left_node.parent = nodes[v]
            nodes[v].left = left_node
            right_node = nodes[val[1]]
            right_node.parent = nodes[v]
            nodes[v].right = right_node

    global cost, chars, rev_chars
    cost = [[0, 1, 1, 1], [1, 0, 1, 1], [1, 1, 0, 1], [1, 1, 1, 0]]
    chars = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    rev_chars = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    for index in range(len(leaf)):
        small_parsimony(root, index)
        # for k, v in nodes.items():
        #     v.vals = [float('inf'), float('inf'), float('inf'), float('inf')]
    return root

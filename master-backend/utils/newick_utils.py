import random


def to_newick_format(labels, result):
    number_sign = 0
    index = 1
    for i in range(0, len(labels)):
        if labels[i] == ')':
            if result == "":
                result = labels[:i + index] + str(number_sign) + labels[i + 1:]
            else:
                result = result[:i + index] + str(number_sign) + labels[i + 1:]
            index += 1
            number_sign += 1
    return result


def tmatrix_to_newick_format(edges):
    result = ";"
    used_nodes = {}
    for (key, value) in edges.items():
        if not used_nodes[key]:
            used_nodes[key] = 1
            result = ")" + str(key) + result


def to_newick(tree):
    newick = ""
    newick = traverse(tree, newick)
    newick = f"{newick};"
    return newick


# take in consideration same name of label!
def traverse(tree, newick):
    if tree.left and not tree.right:
        newick = f"(,{traverse(tree.left, newick)})'{tree.string}'"
    elif not tree.left and tree.right:
        newick = f"({traverse(tree.right, newick)},)'{tree.string}'"
    elif tree.left and tree.right:
        newick = f"({traverse(tree.right, newick)},{traverse(tree.left, newick)})'{tree.string}'"
    elif not tree.left and not tree.right:
        newick = f"'{tree.string}'"
    else:
        pass
    return newick


def to_newick_numbers(tree):
    newick = ""
    newick = traversePhylogeny(tree, newick)
    newick = f"{newick};"
    newick = remove_obsolete_parts(newick)
    return newick


def traversePhylogeny(tree, newick):
    if tree.left and not tree.right:
        newick = f"(,{traversePhylogeny(tree.left, newick)}){tree.string}"
    elif not tree.left and tree.right:
        newick = f"({traversePhylogeny(tree.right, newick)},){tree.string}"
    elif tree.left and tree.right:
        newick = f"({traversePhylogeny(tree.right, newick)},{traversePhylogeny(tree.left, newick)}){tree.string}"
    elif not tree.left and not tree.right:
        newick = f"{tree.string}"
    else:
        pass
    return newick

def remove_obsolete_parts(newick):
    if newick[1] == ',':
        return newick[:1] + newick[2:]
    return newick

def to_newick_parsimony(tree):
    newick = ""
    newick = traverse_parsimony(tree, newick, 0)
    newick = f"{newick};"
    return fix_label_keys(newick)

def fix_label_keys(newick):
    num = 0
    resultString = ""
    for i in range(len(newick)):
        if newick[i].isalnum():
            resultString += newick[i] + str(num)
            num += 1
        else:
            resultString += newick[i]
    return resultString

# take in consideration same name of label!
def traverse_parsimony(tree, newick, index):
    if tree.left and not tree.right:
        newick = f"(,{traverse_parsimony(tree.left, newick, index + 1)}){tree.string}"
    elif not tree.left and tree.right:
        newick = f"({traverse_parsimony(tree.right, newick, index + 1)},){tree.string}"
    elif tree.left and tree.right:
        newick = f"({traverse_parsimony(tree.right, newick, index + 1)},{traverse_parsimony(tree.left, newick, index + 2)}){tree.string}"
    elif not tree.left and not tree.right:
        newick = f"{tree.string}"
    else:
        pass
    return newick

def traverse_parsimony_steps(tree, steps, index):
    if tree.left and not tree.right:
        traverse_parsimony_steps(tree.left, steps, index + 1)
        steps.append({'step': index, 'char': tree.string, 'vals': tree.vals, 'valsCalculations': tree.valsCalculations})
    elif not tree.left and tree.right:
        traverse_parsimony_steps(tree.right, steps, index + 1)
        steps.append({'step': index, 'char': tree.string, 'vals': tree.vals, 'valsCalculations': tree.valsCalculations})
    elif tree.left and tree.right:
        traverse_parsimony_steps(tree.left, steps, index + 1)
        traverse_parsimony_steps(tree.right, steps, index + 1)
        steps.append({'step': index, 'char': tree.string, 'vals': tree.vals, 'valsCalculations': tree.valsCalculations})
    elif not tree.left and not tree.right:
        steps.append({'step': index, 'char': tree.string, 'vals': tree.vals, 'valsCalculations': tree.valsCalculations})
    else:
        pass
    return steps

def sort_parsimony_steps(steps):
    sorted_steps = []
    max_step = find_max_step(steps)
    index = max_step
    while index >= 0:
        for i in range(len(steps)):
            if steps[i]['step'] == index:
                infArray = ""
                for j in range(len(steps[i]['vals']) - 1):
                    infArray += str(steps[i]['vals'][j]) + ','
                infArray += str(steps[i]['vals'][len(steps[i]['vals']) - 1])
                sorted_steps.append((steps[i]['step'], steps[i]['char'], infArray))
        index -= 1
    return sorted_steps

def sort_parsimony_steps_calculations(steps):
    sorted_steps = []
    max_step = find_max_step(steps)
    index = max_step - 1
    while index >= 0:
        for i in range(len(steps)):
            if steps[i]['step'] == index:
                infArray = ""
                for j in range(len(steps[i]['valsCalculations']) - 1):
                    infArray += str(steps[i]['valsCalculations'][j]) + ','
                infArray += str(steps[i]['valsCalculations'][len(steps[i]['valsCalculations']) - 1])
                sorted_steps.append((steps[i]['step'], steps[i]['char'], infArray))
        index -= 1
    return sorted_steps

def find_max_step(steps):
    max = 0
    for i in range(len(steps)):
        if steps[i]['step'] > max:
            max = steps[i]['step']

    return max

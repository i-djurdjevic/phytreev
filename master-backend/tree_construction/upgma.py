# A Quick Implementation of UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
from utils.newick_utils import to_newick_format
import copy
import itertools
import numbers

from Bio.Align import MultipleSeqAlignment
from Bio.Align import substitution_matrices
from Bio.Phylo import BaseTree

class _Matrix:
    def __init__(self, names, matrix=None):

        # check names
        if isinstance(names, list) and all(isinstance(s, str) for s in names):
            if len(set(names)) == len(names):
                self.names = names
            else:
                raise ValueError("Duplicate names found")
        else:
            raise TypeError("'names' should be a list of strings")

        # check matrix
        if matrix is None:
            # create a new one with 0 if matrix is not assigned
            matrix = [[0] * i for i in range(1, len(self) + 1)]
            self.matrix = matrix
        else:
            # check if all elements are numbers
            if (
                    isinstance(matrix, list)
                    and all(isinstance(l, list) for l in matrix)
                    and all(
                isinstance(n, numbers.Number)
                for n in [item for sublist in matrix for item in sublist]
            )
            ):
                # check if the same length with names
                if len(matrix) == len(names):
                    # check if is lower triangle format
                    if [len(m) for m in matrix] == list(range(1, len(self) + 1)):
                        self.matrix = matrix
                    else:
                        raise ValueError("'matrix' should be in lower triangle format")
                else:
                    raise ValueError("'names' and 'matrix' should be the same size")
            else:
                raise TypeError("'matrix' should be a list of numerical lists")

    def __getitem__(self, item):
        # Handle single indexing
        if isinstance(item, (int, str)):
            index = None
            if isinstance(item, int):
                index = item
            elif isinstance(item, str):
                if item in self.names:
                    index = self.names.index(item)
                else:
                    raise ValueError("Item not found.")
            else:
                raise TypeError("Invalid index type.")
            # check index
            if index > len(self) - 1:
                raise IndexError("Index out of range.")
            return [self.matrix[index][i] for i in range(0, index)] + [
                self.matrix[i][index] for i in range(index, len(self))
            ]
        # Handle double indexing
        elif len(item) == 2:
            row_index = None
            col_index = None
            if all(isinstance(i, int) for i in item):
                row_index, col_index = item
            elif all(isinstance(i, str) for i in item):
                row_name, col_name = item
                if row_name in self.names and col_name in self.names:
                    row_index = self.names.index(row_name)
                    col_index = self.names.index(col_name)
                else:
                    raise ValueError("Item not found.")
            else:
                raise TypeError("Invalid index type.")
            # check index
            if row_index > len(self) - 1 or col_index > len(self) - 1:
                raise IndexError("Index out of range.")
            if row_index > col_index:
                return self.matrix[row_index][col_index]
            else:
                return self.matrix[col_index][row_index]
        else:
            raise TypeError("Invalid index type.")

    def __setitem__(self, item, value):
        # Handle single indexing
        if isinstance(item, (int, str)):
            index = None
            if isinstance(item, int):
                index = item
            elif isinstance(item, str):
                if item in self.names:
                    index = self.names.index(item)
                else:
                    raise ValueError("Item not found.")
            else:
                raise TypeError("Invalid index type.")
            # check index
            if index > len(self) - 1:
                raise IndexError("Index out of range.")
            # check and assign value
            if isinstance(value, list) and all(
                    isinstance(n, numbers.Number) for n in value
            ):
                if len(value) == len(self):
                    for i in range(0, index):
                        self.matrix[index][i] = value[i]
                    for i in range(index, len(self)):
                        self.matrix[i][index] = value[i]
                else:
                    raise ValueError("Value not the same size.")
            else:
                raise TypeError("Invalid value type.")
        # Handle double indexing
        elif len(item) == 2:
            row_index = None
            col_index = None
            if all(isinstance(i, int) for i in item):
                row_index, col_index = item
            elif all(isinstance(i, str) for i in item):
                row_name, col_name = item
                if row_name in self.names and col_name in self.names:
                    row_index = self.names.index(row_name)
                    col_index = self.names.index(col_name)
                else:
                    raise ValueError("Item not found.")
            else:
                raise TypeError("Invalid index type.")
            # check index
            if row_index > len(self) - 1 or col_index > len(self) - 1:
                raise IndexError("Index out of range.")
            # check and assign value
            if isinstance(value, numbers.Number):
                if row_index > col_index:
                    self.matrix[row_index][col_index] = value
                else:
                    self.matrix[col_index][row_index] = value
            else:
                raise TypeError("Invalid value type.")
        else:
            raise TypeError("Invalid index type.")

    def __delitem__(self, item):
        index = None
        if isinstance(item, int):
            index = item
        elif isinstance(item, str):
            index = self.names.index(item)
        else:
            raise TypeError("Invalid index type.")
        # remove distances related to index
        for i in range(index + 1, len(self)):
            del self.matrix[i][index]
        del self.matrix[index]
        # remove name
        del self.names[index]

    def insert(self, name, value, index=None):
        if isinstance(name, str):
            # insert at the given index or at the end
            if index is None:
                index = len(self)
            if not isinstance(index, int):
                raise TypeError("Invalid index type.")
            # insert name
            self.names.insert(index, name)
            # insert elements of 0, to be assigned
            self.matrix.insert(index, [0] * index)
            for i in range(index, len(self)):
                self.matrix[i].insert(index, 0)
            # assign value
            self[index] = value
        else:
            raise TypeError("Invalid name type.")

    def __len__(self):
        return len(self.names)

    def __repr__(self):
        return self.__class__.__name__ + "(names=%s, matrix=%s)" % tuple(
            map(repr, (self.names, self.matrix))
        )

    def __str__(self):
        matrix_string = "\n".join(
            [
                self.names[i] + "\t" + "\t".join([str(n) for n in self.matrix[i]])
                for i in range(0, len(self))
            ]
        )
        matrix_string = matrix_string + "\n\t" + "\t".join(self.names)
        return matrix_string

class DistanceMatrix(_Matrix):
    def __init__(self, names, matrix=None):
        _Matrix.__init__(self, names, matrix)
        self._set_zero_diagonal()

    def __setitem__(self, item, value):
        _Matrix.__setitem__(self, item, value)
        self._set_zero_diagonal()

    def _set_zero_diagonal(self):
        for i in range(0, len(self)):
            self.matrix[i][i] = 0

    def format_phylip(self, handle):
        handle.write(f"    {len(self.names)}\n")
        # Phylip needs space-separated, vertically aligned columns
        name_width = max(12, max(map(len, self.names)) + 1)
        value_fmts = ("{" + str(x) + ":.4f}" for x in range(1, len(self.matrix) + 1))
        row_fmt = "{0:" + str(name_width) + "s}" + "  ".join(value_fmts) + "\n"
        for i, (name, values) in enumerate(zip(self.names, self.matrix)):
            # Mirror the matrix values across the diagonal
            mirror_values = (self.matrix[j][i] for j in range(i + 1, len(self.matrix)))
            fields = itertools.chain([name], values, mirror_values)
            handle.write(row_fmt.format(*fields))


# Shim for compatibility with Biopython<1.70 (#1304)
_DistanceMatrix = DistanceMatrix


class DistanceCalculator:
    protein_alphabet = set("ABCDEFGHIKLMNPQRSTVWXYZ")

    dna_models = []
    protein_models = []

    # matrices available
    names = substitution_matrices.load()
    for name in names:
        matrix = substitution_matrices.load(name)
        if name == "NUC.4.4":
            # BLAST nucleic acid scoring matrix
            name = "blastn"
        else:
            name = name.lower()
        if protein_alphabet.issubset(set(matrix.alphabet)):
            protein_models.append(name)
        else:
            dna_models.append(name)

    del protein_alphabet
    del name
    del names
    del matrix

    models = ["identity"] + dna_models + protein_models

    def __init__(self, model="identity", skip_letters=None):
        # Shim for backward compatibility (#491)
        if skip_letters:
            self.skip_letters = skip_letters
        elif model == "identity":
            self.skip_letters = ()
        else:
            self.skip_letters = ("-", "*")

        if model == "identity":
            self.scoring_matrix = None
        elif model in self.models:
            if model == "blastn":
                name = "NUC.4.4"
            else:
                name = model.upper()
            self.scoring_matrix = substitution_matrices.load(name)
        else:
            raise ValueError(
                "Model not supported. Available models: " + ", ".join(self.models)
            )

    def _pairwise(self, seq1, seq2):
        score = 0
        max_score = 0
        if self.scoring_matrix is None:
            # Score by character identity, not skipping any special letters
            score = sum(
                l1 == l2
                for l1, l2 in zip(seq1, seq2)
                if l1 not in self.skip_letters and l2 not in self.skip_letters
            )
            max_score = len(seq1)
        else:
            max_score1 = 0
            max_score2 = 0
            for i in range(0, len(seq1)):
                l1 = seq1[i]
                l2 = seq2[i]
                if l1 in self.skip_letters or l2 in self.skip_letters:
                    continue
                try:
                    max_score1 += self.scoring_matrix[l1, l1]
                except IndexError:
                    raise ValueError(
                        "Bad letter '%s' in sequence '%s' at position '%s'"
                        % (l1, seq1.id, i)
                    ) from None
                try:
                    max_score2 += self.scoring_matrix[l2, l2]
                except IndexError:
                    raise ValueError(
                        "Bad letter '%s' in sequence '%s' at position '%s'"
                        % (l2, seq2.id, i)
                    ) from None
                score += self.scoring_matrix[l1, l2]
            # Take the higher score if the matrix is asymmetrical
            max_score = max(max_score1, max_score2)
        if max_score == 0:
            return 1  # max possible scaled distance
        return 1 - (score * 1.0 / max_score)

    def get_distance(self, msa):
        if not isinstance(msa, MultipleSeqAlignment):
            raise TypeError("Must provide a MultipleSeqAlignment object.")

        names = [s.id for s in msa]
        dm = DistanceMatrix(names)
        for seq1, seq2 in itertools.combinations(msa, 2):
            dm[seq1.id, seq2.id] = self._pairwise(seq1, seq2)
        return dm


class TreeConstructor:

    def build_tree(self, msa):
        raise NotImplementedError("Method not implemented!")

class DistanceCalculator:
    protein_alphabet = set("ABCDEFGHIKLMNPQRSTVWXYZ")

    dna_models = []
    protein_models = []

    # matrices available
    names = substitution_matrices.load()
    for name in names:
        matrix = substitution_matrices.load(name)
        if name == "NUC.4.4":
            # BLAST nucleic acid scoring matrix
            name = "blastn"
        else:
            name = name.lower()
        if protein_alphabet.issubset(set(matrix.alphabet)):
            protein_models.append(name)
        else:
            dna_models.append(name)

    del protein_alphabet
    del name
    del names
    del matrix

    models = ["identity"] + dna_models + protein_models

    def __init__(self, model="identity", skip_letters=None):
        # Shim for backward compatibility (#491)
        if skip_letters:
            self.skip_letters = skip_letters
        elif model == "identity":
            self.skip_letters = ()
        else:
            self.skip_letters = ("-", "*")

        if model == "identity":
            self.scoring_matrix = None
        elif model in self.models:
            if model == "blastn":
                name = "NUC.4.4"
            else:
                name = model.upper()
            self.scoring_matrix = substitution_matrices.load(name)
        else:
            raise ValueError(
                "Model not supported. Available models: " + ", ".join(self.models)
            )

    def _pairwise(self, seq1, seq2):
        score = 0
        max_score = 0
        if self.scoring_matrix is None:
            # Score by character identity, not skipping any special letters
            score = sum(
                l1 == l2
                for l1, l2 in zip(seq1, seq2)
                if l1 not in self.skip_letters and l2 not in self.skip_letters
            )
            max_score = len(seq1)
        else:
            max_score1 = 0
            max_score2 = 0
            for i in range(0, len(seq1)):
                l1 = seq1[i]
                l2 = seq2[i]
                if l1 in self.skip_letters or l2 in self.skip_letters:
                    continue
                try:
                    max_score1 += self.scoring_matrix[l1, l1]
                except IndexError:
                    raise ValueError(
                        "Bad letter '%s' in sequence '%s' at position '%s'"
                        % (l1, seq1.id, i)
                    ) from None
                try:
                    max_score2 += self.scoring_matrix[l2, l2]
                except IndexError:
                    raise ValueError(
                        "Bad letter '%s' in sequence '%s' at position '%s'"
                        % (l2, seq2.id, i)
                    ) from None
                score += self.scoring_matrix[l1, l2]
            # Take the higher score if the matrix is asymmetrical
            max_score = max(max_score1, max_score2)
        if max_score == 0:
            return 1  # max possible scaled distance
        return 1 - (score * 1.0 / max_score)

    def get_distance(self, msa):
        if not isinstance(msa, MultipleSeqAlignment):
            raise TypeError("Must provide a MultipleSeqAlignment object.")

        names = [s.id for s in msa]
        dm = DistanceMatrix(names)
        for seq1, seq2 in itertools.combinations(msa, 2):
            dm[seq1.id, seq2.id] = self._pairwise(seq1, seq2)
        return dm


class UpgmaDistanceTreeConstructor(TreeConstructor):
    methods = ["nj", "upgma"]

    def __init__(self, distance_calculator=None, method="nj"):
        if distance_calculator is None or isinstance(
                distance_calculator, DistanceCalculator
        ):
            self.distance_calculator = distance_calculator
        else:
            raise TypeError("Must provide a DistanceCalculator object.")
        if isinstance(method, str) and method in self.methods:
            self.method = method
        else:
            raise TypeError(
                "Bad method: "
                + method
                + ". Available methods: "
                + ", ".join(self.methods)
            )

    def build_tree(self, msa):
        if self.distance_calculator:
            dm = self.distance_calculator.get_distance(msa)
            tree = None
            if self.method == "upgma":
                tree = self.upgma(dm)
            else:
                tree = self.nj(dm)
            return tree
        else:
            raise TypeError("Must provide a DistanceCalculator object.")

    def upgma(self, distance_matrix):
        # make a copy of the distance matrix to be used
        dm = copy.deepcopy(distance_matrix)
        matrices = []
        matrices.append(distance_matrix)
        branch_lengths = []
        branch_length_calculations = []
        step_distance_calculation = []
        new_distance_calculations = []
        upgma_touples = []
        # init terminal clades
        clades = [BaseTree.Clade(None, name) for name in dm.names]
        # init minimum index
        min_i = 0
        min_j = 0
        inner_count = 0
        while len(dm) > 1:
            min_dist = dm[1, 0]
            # find minimum index
            for i in range(1, len(dm)):
                for j in range(0, i):
                    if min_dist >= dm[i, j]:
                        min_dist = dm[i, j]
                        min_i = i
                        min_j = j

            # create clade
            clade1 = clades[min_i]
            clade2 = clades[min_j]
            inner_count += 1
            inner_clade = BaseTree.Clade(None, "(" + str(clade1.name) + "," + str(clade2.name) + ")")
            inner_clade.clades.append(clade1)
            inner_clade.clades.append(clade2)
            upgma_touples.append((min_dist, clade1.name, clade2.name))

            # assign branch length
            if clade1.is_terminal():
                clade1.branch_length = min_dist * 1.0 / 2
                clade1.branch_length = round(clade1.branch_length, 2)
                clade1_calculation = "d(" + clade1.name + ") / 2 = " + str(min_dist) + "/2 = " + str(clade1.branch_length)
            else:
                clade1.branch_length = min_dist * 1.0 / 2 - self._height_of(clade1)
                clade1.branch_length = round(clade1.branch_length, 2)
                clade1_calculation = "d(" + clade1.name + ") / 2 = " + str(min_dist) + "/2 - " + str(self._height_of(clade1)) + " = " + str(clade1.branch_length)

            if clade2.is_terminal():
                clade2.branch_length = min_dist * 1.0 / 2
                clade2.branch_length = round(clade2.branch_length, 2)
                clade2_calculation = "d(" + clade2.name + ") / 2 = " + str(min_dist) + "/2 = " + str(clade2.branch_length)
            else:
                clade2.branch_length = min_dist * 1.0 / 2 - self._height_of(clade2)
                clade2.branch_length = round(clade2.branch_length, 2)
                clade2_calculation = "d(" + clade2.name + ") / 2 = " + str(min_dist) + "/2 - " + str(self._height_of(clade2)) + " = " + str(
                    clade2.branch_length)

            branch_length_calculations.append((clade1_calculation, clade2_calculation))
            branch_lengths.append((clade1.branch_length, clade1.name, clade2.branch_length, clade2.name))

            # update node list
            clades[min_j] = inner_clade
            del clades[min_i]

            # rebuild distance matrix,
            # set the distances of new node at the index of min_j
            for k in range(0, len(dm)):
                if k != min_i and k != min_j:
                    dm[min_j, k] = (dm[min_i, k] + dm[min_j, k]) * 1.0 / 2
                    dm[min_j, k] = round(dm[min_j, k], 2)
                    distance_calculation = "(d(" + str(dm.names[min_i]) + "," + str(dm.names[k]) + ") + d(" + str(
                        dm.names[min_j]) + "," + str(dm.names[k]) + ")) / 2 = \n" + str(dm[min_i, k]) + " + " + str(dm[min_j, k]) + " / 2 = " + str(dm[min_j, k])
                    step_distance_calculation.append(distance_calculation)

            if step_distance_calculation != []:
                new_distance_calculations.append(step_distance_calculation)
            step_distance_calculation = []

            dm.names[min_j] = "(" + clade1.name + "," + clade2.name + ")"

            del dm[min_i]
            stepMatrix = None
            stepMatrix = copy.deepcopy(dm)
            matrices.append(stepMatrix)
        inner_clade.branch_length = 0
        return (matrices, upgma_touples, BaseTree.Tree(inner_clade), branch_lengths, branch_length_calculations,
                new_distance_calculations)

    def _height_of(self, clade):
        height = 0
        if clade.is_terminal():
            height = clade.branch_length
        else:
            height = height + max(self._height_of(c) for c in clade.clades)
        return height

# lowest_cell:
#   Locates the smallest cell in the table
def lowest_cell(table):
    # Set default to infinity
    min_cell = float("inf")
    x, y = -1, -1

    # Go through every cell, looking for the lowest
    for i in range(len(table)):
        for j in range(len(table[i])):
            if table[i][j] < min_cell:
                min_cell = table[i][j]
                x, y = i, j

    # Return the x, y co-ordinate of cell
    return x, y


# join_labels:
#   Combines two labels in a list of labels
def join_labels(labels, a, b):
    # Swap if the indices are not ordered
    if b < a:
        a, b = b, a

    # Join the labels in the first index
    labels[a] = "(" + labels[a] + "," + labels[b] + ")"

    # Remove the (now redundant) label in the second index
    del labels[b]


# join_table:
#   Joins the entries of a table on the cell (a, b) by averaging their data entries
def join_table(table, a, b):
    # Swap if the indices are not ordered
    if b < a:
        a, b = b, a

    # For the lower index, reconstruct the entire row (A, i), where i < A
    row = []
    for i in range(0, a):
        row.append((table[a][i] + table[b][i]) / 2)
    table[a] = row

    # Then, reconstruct the entire column (i, A), where i > A
    #   Note: Since the matrix is lower triangular, row b only contains values for indices < b
    for i in range(a + 1, b):
        table[i][a] = (table[i][a] + table[b][i]) / 2

    #   We get the rest of the values from row i
    for i in range(b + 1, len(table)):
        table[i][a] = (table[i][a] + table[i][b]) / 2
        # Remove the (now redundant) second index column entry
        del table[i][b]

    # Remove the (now redundant) second index row
    del table[b]


# UPGMA:
#   Runs the UPGMA algorithm on a labelled table
def UPGMA(table, labels):
    # Until all labels have been joined...
    while len(labels) > 1:
        # Locate lowest cell in the table
        x, y = lowest_cell(table)

        # Join the table on the cell co-ordinates
        join_table(table, x, y)

        # Update the labels accordingly
        join_labels(labels, x, y)

    # Return the final label
    return to_newick_format(labels[0], "")

# alpha_labels:
#   Makes labels from a starting letter to an ending letter
def alpha_labels(start, end):
    labels = []
    for i in range(ord(start), ord(end) + 1):
        labels.append(chr(i))
    return labels



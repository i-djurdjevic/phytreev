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


class DistanceTreeConstructor(TreeConstructor):
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
            if self.method == "nj":
                tree = self.nj(dm)
            return tree
        else:
            raise TypeError("Must provide a DistanceCalculator object.")

    def nj(self, distance_matrix):
        nj_touples = []
        matrices = []
        limbLengthArray = []
        # make a copy of the distance matrix to be used
        dm = copy.deepcopy(distance_matrix)
        # matrices = [distance_matrix]
        # init terminal clades
        clades = [BaseTree.Clade(None, name) for name in dm.names]
        # init node distance
        node_dist = [0] * len(dm)
        # init minimum index
        min_i = 0
        min_j = 0
        inner_count = 0
        # special cases for Minimum Alignment Matrices
        if len(dm) == 1:
            root = clades[0]

            return BaseTree.Tree(root, rooted=False)
        elif len(dm) == 2:
            stepMatrix = None
            stepMatrix = copy.deepcopy(dm)
            matrices.append(stepMatrix)
            # minimum distance will always be [1,0]
            min_i = 1
            min_j = 0
            clade1 = clades[min_i]
            clade2 = clades[min_j]
            clade1.branch_length = dm[min_i, min_j] / 2.0
            clade2.branch_length = dm[min_i, min_j] - clade1.branch_length
            inner_clade = BaseTree.Clade(None, "Inner")
            inner_clade.clades.append(clade1)
            inner_clade.clades.append(clade2)
            clades[0] = inner_clade
            root = clades[0]

            return BaseTree.Tree(root, rooted=False)
        while len(dm) > 2:
            stepMatrix = None
            limbLengths = None
            stepMatrix = copy.deepcopy(dm)
            matrices.append(stepMatrix)
            # calculate nodeDist
            for i in range(0, len(dm)):
                node_dist[i] = 0
                for j in range(0, len(dm)):
                    node_dist[i] += dm[i, j]
                node_dist[i] = node_dist[i] / (len(dm) - 2)

            # find minimum distance pair
            min_dist = dm[1, 0] - node_dist[1] - node_dist[0]
            min_i = 0
            min_j = 1
            for i in range(1, len(dm)):
                for j in range(0, i):
                    temp = dm[i, j] - node_dist[i] - node_dist[j]
                    if min_dist > temp:
                        min_dist = temp
                        min_i = i
                        min_j = j
            # create clade
            clade1 = clades[min_i]
            clade2 = clades[min_j]
            inner_count += 1
            inner_clade = BaseTree.Clade(None, "(" + clade1.name + "," + clade2.name + ")")
            inner_clade.clades.append(clade1)
            inner_clade.clades.append(clade2)
            # assign branch length
            clade1.branch_length = (
                                           dm[min_i, min_j] + node_dist[min_i] - node_dist[min_j]
                                   ) / 2.0
            clade2.branch_length = dm[min_i, min_j] - clade1.branch_length

            # we'll append the min distance calculated in the other nj function
            nj_touples.append((clade1.name, clade2.name))

            # nj_touples.append((math.ceil(min_dist), clade1.name, clade2.name))
            # nj_touples.append((clade1.name, clade1.branch_length, clade2.name, clade2.branch_length, min_dist))

            # update node list
            clades[min_j] = inner_clade
            del clades[min_i]

            # rebuild distance matrix,
            # set the distances of new node at the index of min_j
            for k in range(0, len(dm)):
                if k != min_i and k != min_j:
                    dm[min_j, k] = (
                                           dm[min_i, k] + dm[min_j, k] - dm[min_i, min_j]
                                   ) / 2.0

            dm.names[min_j] = inner_clade.name
            del dm[min_i]
        # set the last clade as one of the child of the inner_clade
        root = None
        if clades[0] == inner_clade:
            clades[0].branch_length = 0
            clades[1].branch_length = dm[1, 0]
            clades[0].clades.append(clades[1])
            root = clades[0]
        else:
            clades[0].branch_length = dm[1, 0]
            clades[1].branch_length = 0
            clades[1].clades.append(clades[0])
            root = clades[1]

        return (nj_touples, BaseTree.Tree(root, rooted=False), matrices)

    def _height_of(self, clade):
        height = 0
        if clade.is_terminal():
            height = clade.branch_length
        else:
            height = height + max(self._height_of(c) for c in clade.clades)
        return height


# #################### Tree Scoring and Searching Classes #####################


class Scorer:

    def get_score(self, tree, alignment):
        raise NotImplementedError("Method not implemented!")


class TreeSearcher:
    def search(self, starting_tree, alignment):
        raise NotImplementedError("Method not implemented!")


class NNITreeSearcher(TreeSearcher):
    def __init__(self, scorer):
        if isinstance(scorer, Scorer):
            self.scorer = scorer
        else:
            raise TypeError("Must provide a Scorer object.")

    def search(self, starting_tree, alignment):
        return self._nni(starting_tree, alignment)

    def _nni(self, starting_tree, alignment):
        best_tree = starting_tree
        while True:
            best_score = self.scorer.get_score(best_tree, alignment)
            temp = best_score
            for t in self._get_neighbors(best_tree):
                score = self.scorer.get_score(t, alignment)
                if score < best_score:
                    best_score = score
                    best_tree = t
            # stop if no smaller score exist
            if best_score >= temp:
                break
        return best_tree

    def _get_neighbors(self, tree):
        # make child to parent dict
        parents = {}
        for clade in tree.find_clades():
            if clade != tree.root:
                node_path = tree.get_path(clade)
                # cannot get the parent if the parent is root. Bug?
                if len(node_path) == 1:
                    parents[clade] = tree.root
                else:
                    parents[clade] = node_path[-2]
        neighbors = []
        root_childs = []
        for clade in tree.get_nonterminals(order="level"):
            if clade == tree.root:
                left = clade.clades[0]
                right = clade.clades[1]
                root_childs.append(left)
                root_childs.append(right)
                if not left.is_terminal() and not right.is_terminal():
                    # make changes around the left_left clade
                    # left_left = left.clades[0]
                    left_right = left.clades[1]
                    right_left = right.clades[0]
                    right_right = right.clades[1]
                    # neightbor 1 (left_left + right_right)
                    del left.clades[1]
                    del right.clades[1]
                    left.clades.append(right_right)
                    right.clades.append(left_right)
                    temp_tree = copy.deepcopy(tree)
                    neighbors.append(temp_tree)
                    # neighbor 2 (left_left + right_left)
                    del left.clades[1]
                    del right.clades[0]
                    left.clades.append(right_left)
                    right.clades.append(right_right)
                    temp_tree = copy.deepcopy(tree)
                    neighbors.append(temp_tree)
                    # change back (left_left + left_right)
                    del left.clades[1]
                    del right.clades[0]
                    left.clades.append(left_right)
                    right.clades.insert(0, right_left)
            elif clade in root_childs:
                # skip root child
                continue
            else:
                # method for other clades
                # make changes around the parent clade
                left = clade.clades[0]
                right = clade.clades[1]
                parent = parents[clade]
                if clade == parent.clades[0]:
                    sister = parent.clades[1]
                    # neighbor 1 (parent + right)
                    del parent.clades[1]
                    del clade.clades[1]
                    parent.clades.append(right)
                    clade.clades.append(sister)
                    temp_tree = copy.deepcopy(tree)
                    neighbors.append(temp_tree)
                    # neighbor 2 (parent + left)
                    del parent.clades[1]
                    del clade.clades[0]
                    parent.clades.append(left)
                    clade.clades.append(right)
                    temp_tree = copy.deepcopy(tree)
                    neighbors.append(temp_tree)
                    # change back (parent + sister)
                    del parent.clades[1]
                    del clade.clades[0]
                    parent.clades.append(sister)
                    clade.clades.insert(0, left)
                else:
                    sister = parent.clades[0]
                    # neighbor 1 (parent + right)
                    del parent.clades[0]
                    del clade.clades[1]
                    parent.clades.insert(0, right)
                    clade.clades.append(sister)
                    temp_tree = copy.deepcopy(tree)
                    neighbors.append(temp_tree)
                    # neighbor 2 (parent + left)
                    del parent.clades[0]
                    del clade.clades[0]
                    parent.clades.insert(0, left)
                    clade.clades.append(right)
                    temp_tree = copy.deepcopy(tree)
                    neighbors.append(temp_tree)
                    # change back (parent + sister)
                    del parent.clades[0]
                    del clade.clades[0]
                    parent.clades.insert(0, sister)
                    clade.clades.insert(0, left)
        return neighbors


# ######################## Parsimony Classes ##########################


class ParsimonyScorer(Scorer):
    def __init__(self, matrix=None):
        if not matrix or isinstance(matrix, _Matrix):
            self.matrix = matrix
        else:
            raise TypeError("Must provide a _Matrix object.")

    def get_score(self, tree, alignment):
        # make sure the tree is rooted and bifurcating
        if not tree.is_bifurcating():
            raise ValueError("The tree provided should be bifurcating.")
        if not tree.rooted:
            tree.root_at_midpoint()
        # sort tree terminals and alignment
        terms = tree.get_terminals()
        terms.sort(key=lambda term: term.name)
        alignment.sort()
        if not all(t.name == a.id for t, a in zip(terms, alignment)):
            raise ValueError(
                "Taxon names of the input tree should be the same with the alignment."
            )
        # term_align = dict(zip(terms, alignment))
        score = 0
        for i in range(len(alignment[0])):
            # parsimony score for column_i
            score_i = 0
            # get column
            column_i = alignment[:, i]
            # skip non-informative column
            if column_i == len(column_i) * column_i[0]:
                continue

            # start calculating score_i using the tree and column_i

            # Fitch algorithm without the penalty matrix
            if not self.matrix:
                # init by mapping terminal clades and states in column_i
                clade_states = dict(zip(terms, [{c} for c in column_i]))
                for clade in tree.get_nonterminals(order="postorder"):
                    clade_childs = clade.clades
                    left_state = clade_states[clade_childs[0]]
                    right_state = clade_states[clade_childs[1]]
                    state = left_state & right_state
                    if not state:
                        state = left_state | right_state
                        score_i = score_i + 1
                    clade_states[clade] = state
            # Sankoff algorithm with the penalty matrix
            else:
                inf = float("inf")
                # init score arrays for terminal clades
                alphabet = self.matrix.names
                length = len(alphabet)
                clade_scores = {}
                for j in range(len(column_i)):
                    array = [inf] * length
                    index = alphabet.index(column_i[j])
                    array[index] = 0
                    clade_scores[terms[j]] = array
                # bottom up calculation
                for clade in tree.get_nonterminals(order="postorder"):
                    clade_childs = clade.clades
                    left_score = clade_scores[clade_childs[0]]
                    right_score = clade_scores[clade_childs[1]]
                    array = []
                    for m in range(length):
                        min_l = inf
                        min_r = inf
                        for n in range(length):
                            sl = self.matrix[alphabet[m], alphabet[n]] + left_score[n]
                            sr = self.matrix[alphabet[m], alphabet[n]] + right_score[n]
                            if min_l > sl:
                                min_l = sl
                            if min_r > sr:
                                min_r = sr
                        array.append(min_l + min_r)
                    clade_scores[clade] = array
                # minimum from root score
                score_i = min(array)
                # TODO: resolve internal states
            score = score + score_i
        return score


class ParsimonyTreeConstructor(TreeConstructor):
    def __init__(self, searcher, starting_tree=None):
        self.searcher = searcher
        self.starting_tree = starting_tree

    def build_tree(self, alignment):
        # if starting_tree is none,
        # create a upgma tree with 'identity' scoring matrix
        if self.starting_tree is None:
            dtc = DistanceTreeConstructor(DistanceCalculator("identity"), "upgma")
            self.starting_tree = dtc.build_tree(alignment)
        return self.searcher.search(self.starting_tree, alignment)

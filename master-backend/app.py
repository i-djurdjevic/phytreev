import json
import os

import numpy as np
from flask import Flask, request, jsonify, make_response
from flask_cors import CORS

import tree_construction.small_parsimony as sp
import utils.file_utils as fu
import utils.matrix_utils as mx
import utils.newick_utils as nf
import utils.request_convert_utils as rcu
from tree_construction import additive_phylogeny
from tree_construction import neighbor_joining
from tree_construction.upgma import UpgmaDistanceTreeConstructor
from utils.tree_construction import DistanceTreeConstructor, DistanceMatrix

app = Flask(__name__)
CORS(app)

global filepath
global hasImportedFile
global file_labels
global file_data

hasImportedFile = False
file_labels = []
file_data = []
filepath = ""


@app.route("/nj", methods=['POST'])
def nj():
    global filepath
    global hasImportedFile
    # reset values:
    neighbor_joining.dStarMatrices = []
    neighbor_joining.dMatrices = []
    neighbor_joining.limbLengthTouples = []
    neighbor_joining.minElements = []
    neighbor_joining.stepDistanceCalculations = []
    neighbor_joining.distanceLabels = []

    m, dim, m_labels = rcu.getMatrixAndLabels(request, hasImportedFile, file_labels, file_data)
    distance_matrix = DistanceMatrix(m_labels, m)
    constructor = DistanceTreeConstructor()

    # get phylogenetic tree
    steps, nj_tree, matrices = constructor.nj(distance_matrix)
    matrixStringArray = mx.getMatricesStringArray(matrices)

    # calculate steps
    mm, dimdim = mx.parse_matrix_to_collection(m, dim)
    neighbor_joining.neighborJoining(mm, dimdim, m_labels)

    hasImportedFile = False
    return make_response(jsonify(
        {"steps": neighbor_joining.parse_steps(steps), "tree": nj_tree.format('newick'), "matrices": matrixStringArray,
         "dStarMatrices": neighbor_joining.dStarMatrices, "limbLengths": neighbor_joining.limbLengthTouples,
         "dStarLabels": neighbor_joining.distanceLabels,
         "dStarCalculations": neighbor_joining.stepDistanceCalculations}), 200)


@app.route('/small_parsimony', methods=['POST'])
def small_parsimony():
    global hasImportedFile
    data = request.json
    if hasImportedFile:
        starting_tree = file_data
    else:
        starting_tree = data["data"]
    one_char_tree = mx.parse_starting_tree(starting_tree)
    # full_result_tree = sp.small_parsimony_wrapper(starting_tree)
    one_char_result_tree = sp.small_parsimony_wrapper(one_char_tree)
    steps = nf.traverse_parsimony_steps(one_char_result_tree, [], 0)
    final_steps = nf.sort_parsimony_steps(steps)
    final_steps_calculations = nf.sort_parsimony_steps_calculations(steps)
    one_char_newick_tree = nf.to_newick_parsimony(one_char_result_tree)
    # full_newick_tree = nf.to_newick_parsimony(full_result_tree)
    hasImportedFile = False
    return make_response(
        jsonify({"score": one_char_result_tree.min_par, "tree": one_char_newick_tree, "steps": final_steps, "steps_calculations": final_steps_calculations}), 200)
    # return make_response(jsonify({"score": one_char_result_tree.min_par, "tree": one_char_newick_tree, "full_tree": full_newick_tree, "full_tree_score": full_result_tree.min_par }), 200)


@app.route('/additive_phylogeny', methods=['POST'])
def additive_phylogeny_method():
    global filepath
    global hasImportedFile
    # reset values
    additive_phylogeny.limbLengths = []
    additive_phylogeny.dTrims = []
    additive_phylogeny.dBalds = []
    additive_phylogeny.nodes = []
    additive_phylogeny.limbLengthCaluclations = []
    m, dim, m_labels = rcu.getMatrixAndLabelsPhylogeny(request, hasImportedFile, file_labels, file_data)
    # labels = copy.deepcopy(m_labels)
    edge, weight, _ = additive_phylogeny.additivePhylogeny(m, dim, dim, m_labels)
    parsed_tree = additive_phylogeny.depth_first_search(edge, 0)
    tree_data = additive_phylogeny.toGraph(parsed_tree)
    newick_tree = nf.to_newick_numbers(tree_data)
    # char_tree_data = mx.set_char_data(newick_tree, labels)
    hasImportedFile = False
    return make_response(json.dumps(
        {"tree": newick_tree, "dTrim": additive_phylogeny.dTrims, "dBalds": additive_phylogeny.dBalds,
         "nodes": additive_phylogeny.nodes,
         "limbLengthCalculations": additive_phylogeny.limbLengthCaluclations,
         "limbLengths": additive_phylogeny.limbLengths}, cls=NumpyEncoder), 200)


@app.route('/upgma', methods=['POST'])
def upgma():
    global filepath
    global hasImportedFile
    m, dim, m_labels = rcu.getMatrixAndLabels(request, hasImportedFile, file_labels, file_data)
    constructor = UpgmaDistanceTreeConstructor()
    distance_matrix = DistanceMatrix(m_labels, m)
    matrices, steps, upgma_tree, branch_lengths, branch_length_calculations, new_distance_calculations = constructor.upgma(
        distance_matrix)
    matrixStringArray = mx.getMatricesStringArray(matrices)
    hasImportedFile = False
    return make_response(json.dumps(
        {"steps": steps, "matrix": upgma_tree.format('newick'), "matrices": matrixStringArray, "labels": m_labels,
         "branch_lengths": branch_lengths, "branch_length_calculations": branch_length_calculations,
         "distance_calculations": new_distance_calculations}, default=default_json), 200)


@app.route('/upload_file', methods=['POST'])
def uploadFile():
    global filepath
    global hasImportedFile
    global file_labels
    global file_data
    if request.method == 'POST':
        if request.files:
            uploaded_file = request.files['File']
            type = request.form["type"]
            # change save path accordingly
            filepath = os.path.join("/Users/isidora/Desktop", uploaded_file.filename)
            uploaded_file.save(filepath)
            hasImportedFile = True
            file_labels, file_data = fu.readFile(filepath, type)
    return make_response(jsonify({"matrix": file_data}), 200)


class NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """

    def default(self, obj):
        if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
                            np.int16, np.int32, np.int64, np.uint8,
                            np.uint16, np.uint32, np.uint64)):
            return int(obj)
        elif isinstance(obj, (np.float_, np.float16, np.float32,
                              np.float64)):
            return float(obj)
        elif isinstance(obj, (np.ndarray,)):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def default_json(t):
    return f'{t}'


if __name__ == '__main__':
    app.run()

app.config['FILE_UPLOADS'] = "/Users/isidora/Desktop"

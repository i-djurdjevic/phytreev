import utils.matrix_utils as mx
from tree_construction.upgma import alpha_labels
from Bio import AlignIO
from io import StringIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
import numpy as np

def getMatrixAndLabels(request, hasImportedFile, file_labels, file_data):
    requestData = request.json
    if not hasImportedFile:
        if requestData["type"] == "distance":
            m, dim = mx.parse_matrix(requestData["data"])
            m_labels = alpha_labels("A", chr(ord('A') + dim - 1))
        elif requestData["type"] == "alignment":
            m, dim, m_labels = get_data_from_alignment(requestData["data"])
    else:
        if requestData["type"] == "alignment":
            m, dim, m_labels = get_data_from_alignment(file_data)
        elif requestData["type"] == "distance":
            m, dim = mx.parse_matrix(file_data)
            m_labels = file_labels
    return m, dim, m_labels


def get_data_from_alignment(aln_data):
    handle = StringIO(aln_data)
    aln = AlignIO.read(handle, 'phylip')
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(aln)
    dim = len(distance_matrix)
    # round up values if necessary:
    m = roundUpValuesOfMatrix(distance_matrix.matrix)
    m_labels = distance_matrix.names
    return m, dim, m_labels

def get_data_from_alignment_phylogeny(aln_data):
    handle = StringIO(aln_data)
    aln = AlignIO.read(handle, 'phylip')
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(aln)
    dim = len(distance_matrix)
    # round up values if necessary:
    parsedMatrix = roundUpValuesOfMatrix(distance_matrix.matrix)
    m_labels = distance_matrix.names
    m = fill_up_matrix(parsedMatrix, dim)
    return m, dim, m_labels


def fill_up_matrix(matrix, dim):
    result_matrix = np.zeros((dim, dim))
    for i in range(dim):
        tmp_array = np.zeros(dim)
        for j in range(dim):
            if j <= i:
                tmp_array[j] += round(matrix[i][j], 2)
            else:
                tmp_array[j] += round(matrix[j][i], 2)
        result_matrix[i] += tmp_array
    return result_matrix


def getMatrixAndLabelsPhylogeny(request, hasImportedFile, file_labels, file_data):
    requestData = request.json
    if not hasImportedFile:
        if requestData["type"] == "distance":
            dim, m = mx.process_matrix_from_request(requestData["data"])
            m_labels = alpha_labels("A", chr(ord('A') + dim - 1))
        elif requestData["type"] == "alignment":
            m, dim, m_labels = get_data_from_alignment_phylogeny(requestData["data"])
    else:
        if requestData["type"] == "alignment":
            m, dim, m_labels = get_data_from_alignment_phylogeny(file_data)
        elif requestData["type"] == "distance":
            m_labels = file_labels
            dim, m = mx.process_matrix_from_request(file_data)
    return m, dim, m_labels

def roundUpValuesOfMatrix(matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            matrix[i][j] = round(matrix[i][j], 2)
    return matrix

def getAlignmentAndTree(request, hasImportedData, fileData):
    # if hasImportedData:
    #     return fileData,
    data = request.json
    return data['alignment'], data['tree']
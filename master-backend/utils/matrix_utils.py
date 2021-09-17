import numpy as np

def parse_matrix(data):
  rows = data.split('\n')
  matrix = []
  for row in rows:
    if row != '':
      row_arr = row.split(' ')
      curr_list = []
      for n in row_arr:
        if n == '0':
          curr_list.append(int(n))
          matrix.append(curr_list)
          curr_list = []
        elif n.isdigit():
          curr_list.append(int(n))
  dim = len(matrix)
  return matrix, dim

def parse_string_matrix(data):
  rows = data.split('\n')
  dim = len(rows)
  matrix = []
  for row in rows:
    row_arr = row.split(' ')
    curr_list = []
    for n in row_arr:
      if n.isdigit():
        curr_list.append(int(n))
    if len(curr_list) > 0:
      matrix.append(curr_list)
  return matrix, dim

def process_matrix_from_request(data):
  lines = data.strip('\n').split('\n')
  dim = len(lines)
  matrix = np.zeros((dim, dim), dtype=np.int32)
  i = 0
  for line in lines:
    rowData = np.array(line.split()).astype(np.int32)
    matrix[i] += rowData
    i += 1
  return (dim, matrix)

def parse_request_matrix_to_collection(matrix):
  data = matrix["data"]
  rows = data.split('\n')
  dim = len(rows)
  i = 0
  matrix = {}
  for row in rows:
    row = row.strip()
    row_arr = row.split(' ')
    matrix[i] = {}
    j = 0
    if len(row_arr) > 1:
      for el in row_arr:
        matrix[i][j] = float(el)
        j += 1
      i += 1
  return matrix, dim

def parse_matrix_to_collection(matrix, dim):
  resultMatrix = {}
  for i in range(dim):
    resultMatrix[i] = {}
    flag = False
    for j in range(dim):
      if flag:
        resultMatrix[i][j] = matrix[j][i]
      else:
        resultMatrix[i][j] = matrix[i][j]
      if i == j:
        flag = True
  return resultMatrix, dim

def getMatricesStringArray(matrices):
  matrixStringArray = []
  matrixrez = ""
  for matrix in matrices:
    names = matrix.names
    for name in names:
      matrixrez += " " + name
    matrixrez += "\n"
    namesIndex = 0
    for element in matrix:
      matrixrez += names[namesIndex] + " "
      for i in range(0, len(element)):
        matrixrez += str(element[i]) + " "
      namesIndex += 1
      matrixrez += "\n"
    matrixStringArray.append(matrixrez)
    matrixrez = " "
  return matrixStringArray

def parse_starting_tree(tree):
  result = ""
  tmp_tree = tree.split("\n")
  for i in range(len(tmp_tree) - 1):
    result += tmp_tree[i][:4] + "\n"
  result += tmp_tree[len(tmp_tree) - 1][:4]
  # remove leftover new lines
  return result.strip("\n")

def set_char_data(newick_tree, labels):
  resultNewickString = ""
  charMap = {}
  for i in range(len(labels)):
    charMap[i] = labels[i]
  for i in range(len(newick_tree)):
    if (newick_tree[i].isdigit()):
      resultNewickString += charMap[int(newick_tree[i])]
    else:
      resultNewickString += newick_tree[i]
  return resultNewickString

# def get_matrix_and_labels_from_request(request):
#   global hasImportedFile
#   global filepath
#   if not hasImportedFile:
#     m, dim = parse_matrix(request.json)
#     m_labels = alpha_labels("A", chr(ord('A') + dim - 1))
#   else:
#     m_labels, data = readFile()
#     m, dim = parse_matrix(data)
#     hasImportedFile = False
#     filepath = None
#   return m, dim, m_labels
#
#
# def get_matrices_as_array_of_strings(matrices):
#   matrixStringArray = []
#   # matrixNewickFormat = []
#   matrixrez = ""
#   for matrix in matrices:
#     # matrixNewickFormat.append(([BaseTree.Clade(None, name) for name in matrix.names]).format('newick'))
#     names = matrix.names
#     for name in names:
#       matrixrez += " " + name
#     matrixrez += "\n"
#     namesIndex = 0
#     for element in matrix:
#       matrixrez += names[namesIndex] + " "
#       for i in range(0, len(element)):
#         matrixrez += str(element[i]) + " "
#       namesIndex += 1
#       matrixrez += "\n"
#     matrixStringArray.append(matrixrez)
#     matrixrez = " "
#   return matrixStringArray


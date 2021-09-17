import numpy as np
import pandas as pd
from tree_construction.upgma import alpha_labels

# filepath = None
# hasImportedFile = False

def process_file(filename):
  """
  data format:
  4
  0   13  21  22
  13  0   12  13
  21  12  0   13
  22  13  13  0
  """
  with open(filename) as f:
    for i, line in enumerate(f):
      line = line.strip()
      if i == 0:
        dim = int(line)
        matrix = np.zeros((dim, dim), dtype=np.int32)
      else:
        rowData = np.array(line.split()).astype(np.int32)
        matrix[i - 1] += rowData
  return (dim, matrix)

#  fix read file for phy files (multiple alignments) or txt
def readFile(filepath, type):
  data = ""
  file_extension = filepath.split(".", 1)[1]
  with open(filepath) as file:
    if file_extension == "csv":
      csv_file = pd.read_csv(file, ',')
      m_labels = list(csv_file.columns)
      for index, row in csv_file.iterrows():
        for label in m_labels:
          data += str(row[label]) + " "
        data += "\n"
    elif file_extension == "txt" and not type == "parsimony":
      lines = file.readlines()
      m_labels = lines[0].strip().split(" ")
      for i in range(1, len(lines)):
        data += lines[i]
    elif file_extension == "txt" and type == "parsimony":
      lines = file.readlines()
      for i in range(len(lines)):
        data += lines[i]
      m_labels = None
    elif file_extension == "phy":
      lines = file.readlines()
      for line in lines:
        data += line
      dim = data.split(" ")[0]
      m_labels = alpha_labels("A", chr(ord('A') + int(dim) - 1))
  return m_labels, data

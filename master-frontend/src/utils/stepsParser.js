import React from "react";

export function parseSteps(data) {
  let result = "";
  for (let i = 0; i < data.length; i++) {
    result += i + 1 + ") ";
    result += data[i];
    result += "\n";
  }
  return result;
}

export function parseSmallParsimonySteps(result, steps, stepsCalculations, step) {
  let resultStepString = "Best parsimony score is: " + result.score + "\n\n";
  if (parseInt(step) === 0) {
    resultStepString += "Initially we set inf in leaves where the char != [A, C, G, T], and 0 where it is:\n\n"
  }
  for (let i = 0; i < steps[step].length; i++) {
    resultStepString += steps[step][i][1] + " -> " + steps[step][i][2] + "\n"
  }
  resultStepString += "\n";
  if (parseInt(step) !== 0) {
    let calcs = stepsCalculations[step - 1];
    for (let j = 0; j < calcs.length; j++) {
      let stepCalcs = calcs[j][2].split(",");
      resultStepString += calcs[j][1] + ":\n";
      for (let k = 0; k < stepCalcs.length; k++) {
        resultStepString += stepCalcs[k] + "\n";
      }
      resultStepString += "\n";
    }
  }
  return resultStepString;
}

export function parseAdditivePhylogenySteps(result, step) {
  let resultStepString = "Smallest calculated limb length is: " + result.limbLengths[step] + "\n\n";
  let limbLengthCalculations = result.limbLengthCalculations[step];
  if (limbLengthCalculations.length > 0) {
    resultStepString += "Limb length calculations were:\n";
    for (let i = 0; i < limbLengthCalculations.length; i++) {
      resultStepString += limbLengthCalculations[i] + "\n";
    }
    if (step == result.limbLengthCalculations.length - 1) {
      resultStepString += "\nRemaining limb length to be calculated: \n" + result.limbLengths[result.limbLengthCalculations.length];
    }
  }
  return resultStepString;
}

export function parseUpgmaSteps(result, step) {
  let smallestElementArray = result.steps[step];
  let resultStepString = "Smallest element in matrix is: "
      + smallestElementArray[0] + " for nodes [" + smallestElementArray[1]
      + "] and " + "[" + smallestElementArray[2] + "]";
  resultStepString += "\n\nBranch lengths:\n"
  for (let i = 0; i < result.branch_length_calculations[step].length; i++) {
    resultStepString += result.branch_length_calculations[step][i] + "\n";
  }
  if (result.distance_calculations.length > step) {
    resultStepString += "\n\nCalculating new elements in the matrix:\n"
    for (let i = 0; i < result.distance_calculations[step].length; i++) {
      resultStepString += result.distance_calculations[step][i] + "\n\n";
    }
  }
  return resultStepString;
}

export function parseNjSteps(result, step) {
  let smallestElementArray = result.steps[step];
  let resultStepString = "Smallest element in matrix is: "
      + smallestElementArray[0] + " for nodes [" + smallestElementArray[1]
      + "] and " + "[" + smallestElementArray[2] + "]\n\n";

  let calculationStrings = result.dStarCalculations[step];
  resultStepString += "D* calculations: \n";
  for (let i = 0; i < calculationStrings.length; i++) {
    resultStepString += calculationStrings[i] + "\n\n";
  }
  return resultStepString;
}

// fix method, remove dStarMatrices as it's redundant (we have that info in result)
export function parseNjMatrix(result, dStarMatrices, step) {
  let stepMatrix = dStarMatrices[step];
  const labels = getLabels(result, step);
  let keys = Object.keys(stepMatrix);
  return <div>
      <p style={{'font-weight':'bold', 'font-size': '15px'}}>D* matrix:</p>
      <table style={{'border':'1px solid black', 'width':'100%'}}>
        <tr style={{'border':'1px solid black'}}>
          {labels.map((item) => <td style={{'border':'1px solid black', 'text-align':'center'}}>{item}</td>)}
        </tr>
      {keys.map((value, index) => {
        return <tr style={{'border':'1px solid black'}}>
            {Object.entries(stepMatrix[value]).map((key, item) => {
              if (item === 0) {
                let vrednosti = [labels[index + 1], key[1]];
                return vrednosti.map(item => <td style={{'border':'1px solid black', 'text-align':'center'}}>{item}</td>)
              } else {
                return <td style={{'border':'1px solid black', 'text-align':'center'}}>{key[1]}</td>
              }
            })}
          </tr>
      }
      )}
    </table>
  </div>;
}

export function parseTable(result, step) {
  let matrixElements = getMatrixElements(result, step);
  return <div>
    <p style={{'font-weight':'bold', 'font-size': '15px'}}>D matrix:</p>
  <table style={{'border':'1px solid black', 'width':'100%'}}>
    {matrixElements.map(value =>
        <tr style={{'border':'1px solid black'}}>
          {value.map((value) => value.map(value => <td style={{'border':'1px solid black', 'text-align':'center'}}>{value}</td>))}
        </tr>
    )}
  </table>
  </div>;
}

export function parseArrayTable(matrix, describeMatrix) {
  return <div>
    <p style={{'font-weight':'bold', 'font-size': '15px'}}>{describeMatrix}</p>
    <table style={{'border':'1px solid black', 'width':'100%'}}>
      {matrix.map(value =>
          <tr style={{'border':'1px solid black'}}>
            {value.map((value) => <td style={{'border':'1px solid black', 'text-align':'center'}}>{value}</td>)}
          </tr>
      )}
    </table>
  </div>;
}

export function parseInitialMatrix(matrixString, result) {
  let matrixArray = []
  let labels = []
  if (result.dBalds) {
    labels = result.dBalds[0][0];
  } else if (result.matrices[0]) {
    labels = getLabels(result, 0);
  } else if (result.labels) {
    labels.push(" ");
    labels.push(result.labels);
  }
  matrixArray.push(labels);
  let matrixRows = matrixString.trim().split("\n");
  for (let i = 0; i < matrixRows.length; i++) {
    let matrixElements = matrixRows[i].trim().split(' ');
    let matrixRow = [];
    matrixRow.push(labels[i + 1]);
    for (let j = 0; j < matrixElements.length; j++) {
      matrixRow.push(matrixElements[j]);
    }
    matrixArray.push(matrixRow);
  }
  return matrixArray;
}

export function parseInitialParsimonySteps(steps) {
  let resultSteps = {}
  let currentStep = steps[0][0];
  let currentStepArray = []
  let index = 0
  for (let i = 0; i < steps.length; i++) {
    if (currentStep == steps[i][0]) {
      currentStepArray.push(steps[i])
      }
    else {
      resultSteps[index] = currentStepArray;
      currentStepArray = [];
      currentStepArray.push(steps[i])
      currentStep = steps[i][0];
      index += 1;
    }
  }
  if (currentStepArray.length !== 0) {
    resultSteps[index] = currentStepArray;
  }
  return resultSteps;
}

const getMatrixElements = (data, step) => {
  let outputMatrix = [];
  let currentArray = [];
  let elements = data.matrices[step].split("\n");
  let labels = elements[0].trim().split(" ");
  labels.unshift(" ");
  currentArray.push(labels);
  outputMatrix.push(currentArray);
  currentArray = [];
  // skipping first row as it only contains labels of the matrix
  for (let i = 1; i < elements.length; i++) {
    let row = elements[i].split(" ").slice(0, -1);
    currentArray.push(row);
    outputMatrix.push(currentArray);
    currentArray = [];
  }
  return outputMatrix;
}

const getLabels = (data, step) => {
  let elements = data.matrices[step].split("\n");
  let labels = elements[0].trim().split(" ");
  labels.unshift(" ");
  return labels;
}

import {
  // postAlignmentAndTree,
  postData
} from "./fetchUtils";

const MAIN_ENDPOINT = "http://localhost:5000/";

export function applyAlgorithm(url, data, type) {
  return postData(MAIN_ENDPOINT + url, data, type);
}


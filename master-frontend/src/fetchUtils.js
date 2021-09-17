function getOptions(data, type) {
  let requestData = {"data": data, "type": type};
  return {
    method: 'POST',
    body: JSON.stringify(requestData),
    headers: {
      'Content-Type': 'application/json'
    }
  };
}

export async function postData(url, data, type) {
  const options = getOptions(data, type);
  const response = await fetch(url, options);
  return await response.json();
}
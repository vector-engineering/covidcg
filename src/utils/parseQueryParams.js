/*
Currently using URLSearchParams()
export function parseQueryParams() {
  var queryString = window.location.search;
  if (queryString === '') return {};
  var params = {};
  queryString = decodeURI(queryString);

  console.log(queryString);

  for (const pairNum in queryString.split('&')) {
    const key = pairNum.split('=')[0];
    if (!key.length) continue;
    if (typeof params[key] === 'undefined') {
      params[key] = [];
    }

    let values = pairNum.split('=')[1];
    values = values.split('%2C');

    if (values.length > 1) {
      if (key === 'customCoordinates' || key === 'residueCoordinates') {
        const oldVals = values;
        const pair = [parseInt(oldVals[0]), parseInt(oldVals[1])];
        values = [];
        values.push(pair);
      }
      params[key] = values;
    } else {
      params[key].push(values);
    }
  }
  return params;
}
*/

export function updateURLFromParams(params) {
  var queryString = '?';
  for (const pair of params.entries()) {
    queryString += pair[0] + '=' + pair[1] + '&';
  }
  queryString = queryString.slice(0, -1);
  window.history.replaceState({}, '', `${location.pathname}${queryString}`);
}

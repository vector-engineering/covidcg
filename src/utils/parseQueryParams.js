export function parseQueryParams() {
  var queryString = window.location.search;
  if (queryString === '') return {};
  var keyValPairs = [];
  var params = {};
  queryString = queryString.replace("?", "");
  queryString = queryString.replace(/%20/g, " ");
  queryString = queryString.replace(/%2C/g, ",")

  if (queryString.length) {
    keyValPairs = queryString.split('&');
    for (var pairNum in keyValPairs) {
      var key = keyValPairs[pairNum].split('=')[0];
      if (!key.length) continue;
      if (typeof params[key] === 'undefined') {
        params[key] = [];
      }

      var values = keyValPairs[pairNum].split('=')[1];
      values = values.split(',');

      if (Array.isArray(values) && values.length > 1) {
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
    }
  return params;
}

export function updateURLFromParams(params) {
  var queryString = '?';
  Object.keys(params).forEach((key) => {
    queryString += key + '=' + params[key] + '&';
  });
  queryString = queryString.slice(0, -1);
  queryString = queryString.replace(/%2C/g, ",")
  window.history.replaceState({}, '', `${location.pathname}${queryString}`)
}

export function updateURLFromParams(params) {
  var queryString = '?';
  for (const pair of params.entries()) {
    queryString += pair[0] + '=' + pair[1] + '&';
  }
  queryString = queryString.slice(0, -1);
  window.history.replaceState({}, '', `${location.pathname}${queryString}`);
}

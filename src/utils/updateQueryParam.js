/* 
export function updateQueryStringParam(key, value) {
  const baseUrl = [
    location.protocol,
    '//',
    location.host,
    location.pathname,
  ].join('');
  const urlQueryString = document.location.search;
  const newParam = key + '=' + value;
  let params = '?' + newParam;

  // If the "search" string exists, then build params from it
  if (urlQueryString) {
    const keyRegex = new RegExp('([?&])' + key + '[^&]*');
    // If param exists already, update it
    if (urlQueryString.match(keyRegex) !== null) {
      params = urlQueryString.replace(keyRegex, '$1' + newParam);
    } else {
      // Otherwise, add it to end of query string
      params = urlQueryString + '&' + newParam;
    }
  }
  window.history.replaceState({}, '', baseUrl + params);
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

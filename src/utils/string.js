// https://flaviocopes.com/how-to-uppercase-first-letter-javascript/
export function capitalize(s) {
  if (typeof s !== 'string') return '';
  return s.charAt(0).toUpperCase() + s.slice(1);
}

// https://stackoverflow.com/questions/6122571/simple-non-secure-hash-function-for-javascript
export function hashCode(str) {
  var hash = 0;
  if (str.length == 0) {
    return hash;
  }
  for (var i = 0; i < str.length; i++) {
    var char = str.charCodeAt(i);
    hash = (hash << 5) - hash + char;
    hash = hash & hash; // Convert to 32bit integer
  }
  return hash;
}

const basePairMap = {
  A: 'T',
  T: 'A',
  C: 'G',
  G: 'C',
};
export function reverseComplement(str) {
  if (typeof str === 'object') {
    const res = {};
    Object.keys(str).forEach((key) => {
      res[key] = Array.from(str[key])
        .reverse()
        .map((char) => basePairMap[char])
        .join('');
    });
    return res;
  } else {
    return Array.from(str)
      .reverse()
      .map((char) => basePairMap[char])
      .join('');
  }
}

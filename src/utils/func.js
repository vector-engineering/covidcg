// Function utilities

// https://github.com/you-dont-need/You-Dont-Need-Lodash-Underscore#_throttle
export function throttle(func, timeFrame) {
  var lastTime = 0;
  return function (...args) {
    var now = new Date();
    if (now - lastTime >= timeFrame) {
      func(...args);
      lastTime = now;
    }
  };
}

// https://github.com/you-dont-need/You-Dont-Need-Lodash-Underscore#_debounce
export function debounce(func, wait, immediate) {
  var timeout;
  return function () {
    var context = this,
      args = arguments;
    clearTimeout(timeout);
    timeout = setTimeout(function () {
      timeout = null;
      if (!immediate) func.apply(context, args);
    }, wait);
    if (immediate && !timeout) func.apply(context, args);
  };
}

// Shallow array equality
export function arrayEqual(a, b) {
  if (a.length !== b.length) return false;
  return a.every((itemA) => b.includes(itemA));
}

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

// Shallow array equality
export function arrayEqual(a, b) {
  if (a.length !== b.length) return false;
  return a.every((itemA) => b.includes(itemA));
}

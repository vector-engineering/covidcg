// NaN-sensitive version of Math.max() for two params
export function nanmax(a, b) {
  if (Number.isNaN(a)) {
    return b;
  } else if (Number.isNaN(b)) {
    return a;
  }
  return Math.max(a, b);
}

// NaN-sensitive version of Math.min() for two params
export function nanmin(a, b) {
  if (Number.isNaN(a)) {
    return b;
  } else if (Number.isNaN(b)) {
    return a;
  }
  return Math.min(a, b);
}

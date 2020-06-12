// Date utilities

export function intToISO(dateNum) {
  return new Date(dateNum).toISOString().substring(0, 10);
}

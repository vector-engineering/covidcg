// Date utilities

export function intToISO(dateNum) {
  return new Date(dateNum).toISOString().substring(0, 10);
}

export function ISOToInt(dateStr) {
  return new Date(dateStr).getTime();
}

export function todayISO() {
  return intToISO(new Date().getTime());
}

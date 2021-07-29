import { DNA_OR_AA, GROUPS } from '../constants/defs.json';

export const formatSnv = (snvStr, dnaOrAa) => {
  // Don't do this if it's a special group
  if (Object.values(GROUPS).includes(snvStr)) {
    return snvStr;
  }

  // Print as REF POS ALT
  // i.e., 23403|A|G -> A23403G, S|614|D|G -> S:D614G
  const chunks = snvStr.split('|');
  if (dnaOrAa === DNA_OR_AA.DNA) {
    return `${chunks[1]}${chunks[0]}${chunks[2]}`;
  } else if (dnaOrAa === DNA_OR_AA.AA) {
    return `${chunks[0]}:${chunks[2]}${chunks[1]}${chunks[3]}`;
  }
};

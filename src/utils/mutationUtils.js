import { DNA_OR_AA, GROUPS } from '../constants/defs.json';

export const formatMutation = (mutationStr, dnaOrAa) => {
  // Don't do this if it's a special group
  if (Object.values(GROUPS).includes(mutationStr)) {
    return mutationStr;
  }

  // Print as REF POS ALT
  // i.e., 23403|A|G -> A23403G, S|614|D|G -> S:D614G

  const chunks = mutationStr.split('|');
  if (dnaOrAa === DNA_OR_AA.DNA) {
    return `${chunks[1]}${chunks[0]}${chunks[2]}`;
  } else if (dnaOrAa === DNA_OR_AA.AA) {
    let ref = chunks[2];
    let alt = chunks[3];

    // Translate stop codon from "_" to more conventional "*"
    if (ref === '_') {
      ref = '*';
    }
    if (alt === '_') {
      alt = '*';
    }

    // For deletions, format as ∆[ref][pos]
    // i.e., F13- --> ∆F13
    if (alt === '-') {
      alt = '';
      ref = '∆' + ref;
    }

    // For insertions, remove ref '-'
    // i.e., -13F --> 13F
    if (ref === '-') {
      ref = '';
    }

    return `${chunks[0]}:${ref}${chunks[1]}${alt}`;
  }
};

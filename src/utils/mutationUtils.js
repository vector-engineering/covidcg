import { config } from '../config';
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
    // If the virus has more than one segment, then append the segment name
    // to the beginning of each mutation
    if (config.segments.length > 1) {
      return `${chunks[0]}:${chunks[2]}${chunks[1]}${chunks[3]}`;
    } else {
      return `${chunks[2]}${chunks[1]}${chunks[3]}`;
    }
  } else if (dnaOrAa === DNA_OR_AA.AA) {
    let feature = chunks[0];
    let pos = parseInt(chunks[1]);
    let ref = chunks[2];
    let alt = chunks[3];

    // For mutations with multiple reference AAs affected, display
    // the position as a range of affected positions
    // e.g., LPPA24S --> LPPA24/27S
    if (ref.length > 1) {
      pos = `${pos}/${pos + ref.length - 1}`;
    }

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

    return `${feature}:${ref}${pos}${alt}`;
  }
};

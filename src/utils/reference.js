import { DEGENERATES } from '../constants/defs.json';
import { reverseComplement } from './string';
// eslint-disable-next-line import/no-unresolved
import references from '../../static_data/__VIRUS__/reference.json';

export function getReferences() {
  return references;
}
export function getReference(referenceName) {
  return references[referenceName];
}

Object.keys(references).forEach((referenceName) => {
  Object.keys(references[referenceName]['segments']).forEach((segment) => {
    references[referenceName]['segments'][segment]['sequence_rc'] =
      reverseComplement(
        references[referenceName]['segments'][segment]['sequence']
      );
  });
});

export function getReferenceNames() {
  return Object.keys(references);
}

const _subtypes = new Set();
const subtypeReferenceMap = {};
Object.keys(references).forEach((referenceName) => {
  const reference = references[referenceName];
  const subtype = reference.subtype;

  // Add subtype to set of all subtypes
  _subtypes.add(subtype);

  // Add subtype-reference mapping
  if (!Object.prototype.hasOwnProperty.call(subtypeReferenceMap, subtype)) {
    subtypeReferenceMap[subtype] = [];
  }
  subtypeReferenceMap[subtype].push(referenceName);
});
const subtypes = Array.from(_subtypes);

export function getSubtypes() {
  return subtypes;
}
export function getSubtypeReferenceMap() {
  return subtypeReferenceMap;
}
export function getReferencesForSubtype(subtype) {
  return subtypeReferenceMap[subtype];
}

export const queryReferenceSequence = (referenceName, query) => {
  let degenInds = [];
  [...query].forEach((char, i) => {
    // Check if char represents a nucleotide (NT)
    if (char.toUpperCase() in DEGENERATES) {
      // If char is a degenerate NT, save an array where:
      // The first element is the index of the degenerate NT in seq
      // The second element will be used to set a valid NT
      degenInds.push([i, 0]);
    } else if (['A', 'C', 'T', 'G'].includes(char.toUpperCase())) {
      // Char won't need to be processed, continue to next char
      return;
    } else {
      // Invalid sequence
      return 0;
    }
  });

  // If no degenerate nucleotides, keep original logic
  if (degenInds.length === 0) {
    return lookForSeq(referenceName, query);
  }

  // If there are degenerate NTs, generate and check the sequences
  // Create the first sequence
  let newseq = query;
  // Calculate the max number of permutations to use later as an e-stop
  let maxPerms = 1;
  degenInds.forEach((item) => {
    // Replace each degenerate with a valid nucleotide
    let validNTArray = DEGENERATES[query.charAt(item[0])];
    newseq = replaceChar(newseq, item[0], validNTArray[item[1]]);
    // Number.MAX_VALUE = 2^1024 so maxPerms shouldn't overflow
    maxPerms *= validNTArray.length;
  });

  // Check if newseq is in the reference
  // (res is an array if newseq is in the reference and 0 otherwise)
  let res = lookForSeq(newseq);
  let counter = 1;
  let degenInd = 0;
  let validNTArray = DEGENERATES[query.charAt(degenInds[degenInd][0])];
  while (res === 0 && counter < maxPerms) {
    // While newseq is not in the reference, generate and check next sequence
    // To generate sequence, change the first degenerate to its next valid NT
    degenInds[degenInd][1]++;
    while (degenInds[degenInd][1] >= validNTArray.length) {
      // If the index is beyond validNTArray reset it to its first valid NT
      degenInds[degenInd][1] = 0;
      newseq = replaceChar(newseq, degenInds[degenInd][0], validNTArray[0]);
      // Change the next degenerate to its next valid NT.
      degenInd++;
      degenInds[degenInd][1]++;
      validNTArray = DEGENERATES[query.charAt(degenInds[degenInd][0])];
      // If this index is beyond validNTArray, repeat.
    }
    // Construct the new sequence
    newseq = replaceChar(
      newseq,
      degenInds[degenInd][0],
      validNTArray[degenInds[degenInd][1]]
    );
    // Increment the permutation counter.
    counter++;
    // Reset to the first degenerate.
    degenInd = 0;
    validNTArray = DEGENERATES[query.charAt(degenInds[degenInd][0])];
    // Check the new seq
    res = lookForSeq(referenceName, newseq);
  }
  return res;
};

function replaceChar(str, ind, char) {
  return str.substr(0, ind) + char + str.substr(ind + 1);
}

// Look for a query sequence in the specified reference
// Returns: array of length 3: [segment, start, end]
const lookForSeq = (referenceName, query) => {
  for (const segment in references[referenceName]['segments']) {
    const sequenceForward =
      references[referenceName]['segments'][segment]['sequence'];
    const sequenceReverse =
      references[referenceName]['segments'][segment]['sequence_rc'];

    const forwardMatch = sequenceForward.indexOf(query);
    const reverseMatch = sequenceReverse.indexOf(query);

    if (forwardMatch > -1) {
      return [segment, forwardMatch + 1, forwardMatch + query.length];
    } else if (reverseMatch > -1) {
      return [
        segment,
        sequenceForward.length - reverseMatch - query.length + 1,
        sequenceForward.length - reverseMatch,
      ];
    }
  }
  return 0;
};

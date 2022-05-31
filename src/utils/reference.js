import { DEGENERATES } from '../constants/defs.json';
import { reverseComplement } from './string';
import { memoize } from './func';
// eslint-disable-next-line import/no-unresolved
import references from '../../static_data/__VIRUS__/reference.json';

export function getReferences() {
  return references;
}
export function getReference(referenceName) {
  return references[referenceName];
}

Object.keys(references).forEach((referenceName) => {
  references[referenceName]['sequence_rc'] = reverseComplement(
    references[referenceName]['sequence']
  );
});

export function getReferenceNames() {
  return Object.keys(references);
}

export function referenceSequenceIncludes(referenceName, query) {
  return (
    forRefSeqIncludes(referenceName, query) ||
    revRefSeqIncludes(referenceName, query)
  );
}

const forRefSeqIncludes = memoize((referenceName, query) => {
  return references[referenceName]['sequence'].includes(query);
});
const revRefSeqIncludes = memoize((referenceName, query) => {
  return references[referenceName]['sequence_rc'].includes(query);
});
const forRefSeqIndexOf = memoize((referenceName, query) => {
  return references[referenceName]['sequence'].indexOf(query);
});
const revRefSeqIndexOf = memoize((referenceName, query) => {
  return references[referenceName]['sequence_rc'].indexOf(query);
});

export const queryReferenceSequence = memoize((referenceName, query) => {
  let ind = [];
  [...query].forEach((char, i) => {
    // Check if char represents a nucleotide (NT)
    if (char.toUpperCase() in DEGENERATES) {
      // If char is a degenerate NT, save an array where:
      // The first element is the index of the degenerate NT in seq
      // The second element will be used to set a valid NT
      ind.push([i, 0]);
    } else if (['A', 'C', 'T', 'G'].includes(char.toUpperCase())) {
      // Char won't need to be processed, continue to next char
      return;
    } else {
      // Invalid sequence
      return 0;
    }
  });

  if (ind.length) {
    // If there are degenerate NTs, generate and check the sequences
    // Create the first sequence
    let newseq = query;
    // Calculate the max number of permutations to use later as an e-stop
    let maxPerms = 1;
    ind.forEach((item) => {
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
    let validNTArray = DEGENERATES[query.charAt(ind[degenInd][0])];
    while (res === 0 && counter < maxPerms) {
      // While newseq is not in the reference, generate and check next sequence
      // To generate sequence, change the first degenerate to its next valid NT
      ind[degenInd][1]++;
      while (ind[degenInd][1] >= validNTArray.length) {
        // If the index is beyond validNTArray reset it to its first valid NT
        ind[degenInd][1] = 0;
        newseq = replaceChar(newseq, ind[degenInd][0], validNTArray[0]);
        // Change the next degenerate to its next valid NT.
        degenInd++;
        ind[degenInd][1]++;
        validNTArray = DEGENERATES[query.charAt(ind[degenInd][0])];
        // If this index is beyond validNTArray, repeat.
      }
      // Construct the new sequence
      newseq = replaceChar(
        newseq,
        ind[degenInd][0],
        validNTArray[ind[degenInd][1]]
      );
      // Increment the permutation counter.
      counter++;
      // Reset to the first degenerate.
      degenInd = 0;
      validNTArray = DEGENERATES[query.charAt(ind[degenInd][0])];
      // Check the new seq
      res = lookForSeq(referenceName, newseq);
    }
    return res;
  } else {
    // No degenerate nucleotides, keep original logic
    return lookForSeq(referenceName, query);
  }
});

function replaceChar(str, ind, char) {
  return str.substr(0, ind) + char + str.substr(ind + 1);
}

const lookForSeq = memoize((referenceName, query) => {
  if (forRefSeqIncludes(referenceName, query)) {
    return [
      forRefSeqIndexOf(referenceName, query) + 1,
      forRefSeqIndexOf(referenceName, query) + query.length,
    ];
  } else if (revRefSeqIncludes(referenceName, query)) {
    return [
      references[referenceName]['sequence'].length -
        revRefSeqIndexOf(referenceName, query) -
        query.length +
        1,
      references[referenceName]['sequence'].length -
        revRefSeqIndexOf(referenceName, query),
    ];
  } else {
    return 0;
  }
});

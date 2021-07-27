import refSeq from '../../static_data/reference.json';
import { DEGENERATES } from '../constants/defs.json';

import { reverseComplement } from './string';
import { memoize } from './func';

const reverseComplementRefSeq = reverseComplement(refSeq.ref_seq);

export function getReferenceSequence() {
  return refSeq.ref_seq;
}
export function getReferenceSequenceReverseComplement() {
  return reverseComplementRefSeq;
}

export function referenceSequenceIncludes(seq) {
  return forRefSeqIncludes(seq) || revRefSeqIncludes(seq);
}

const forRefSeqIncludes = memoize((seq) => {
  return refSeq.ref_seq.includes(seq);
});
const revRefSeqIncludes = memoize((seq) => {
  return reverseComplementRefSeq.includes(seq);
});

export const queryReferenceSequence = memoize((seq) => {
  let ind = [];
  [...seq].forEach((char, i) => {
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
    let newseq = seq;
    // Calculate the max number of permutations to use later as an e-stop
    let maxPerms = 1;
    ind.forEach((item) => {
      // Replace each degenerate with a valid nucleotide
      let validNTArray = DEGENERATES[seq.charAt(item[0])];
      newseq = replaceChar(newseq, item[0], validNTArray[item[1]]);
      // Number.MAX_VALUE = 2^1024 so maxPerms shouldn't overflow
      maxPerms *= validNTArray.length;
    });

    // Check if newseq is in the reference
    // (res is an array if newseq is in the reference and 0 otherwise)
    let res = lookForSeq(newseq);
    let counter = 1;
    let degenInd = 0;
    let validNTArray = DEGENERATES[seq.charAt(ind[degenInd][0])];
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
        validNTArray = DEGENERATES[seq.charAt(ind[degenInd][0])];
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
      validNTArray = DEGENERATES[seq.charAt(ind[degenInd][0])];
      // Check the new seq
      res = lookForSeq(newseq);
    }
    return res;
  } else {
    // No degenerate nucleotides, keep original logic
    return lookForSeq(seq);
  }
});

const forRefSeqIndexOf = memoize((seq) => {
  return refSeq.ref_seq.indexOf(seq);
});
const revRefSeqIndexOf = memoize((seq) => {
  return reverseComplementRefSeq.indexOf(seq);
});

function replaceChar(str, ind, char) {
  return str.substr(0, ind) + char + str.substr(ind + 1);
}

const lookForSeq = memoize((seq) => {
  if (forRefSeqIncludes(seq)) {
    return [forRefSeqIndexOf(seq) + 1, forRefSeqIndexOf(seq) + seq.length];
  } else if (revRefSeqIncludes(seq)) {
    return [
      refSeq.ref_seq.length - revRefSeqIndexOf(seq) - seq.length + 1,
      refSeq.ref_seq.length - revRefSeqIndexOf(seq),
    ];
  } else {
    return 0;
  }
});

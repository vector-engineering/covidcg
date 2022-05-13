import { DEGENERATES } from '../constants/defs.json';
import { reverseComplement } from './string';
import { memoize } from './func';
import refSeq from '../../static_data/__VIRUS__/reference.json';

const reverseComplementRefSeq = Object.assign({}, refSeq);
Object.keys(reverseComplementRefSeq).forEach((key) => {
  reverseComplementRefSeq[key] = reverseComplement(
    reverseComplementRefSeq[key]
  );
});

export function getReferenceSequence(key = '') {
  return key ? refSeq.ref_seq[key] : refSeq.ref_seq;
}
export function getReferenceSequenceReverseComplement(key = '') {
  return key ? reverseComplementRefSeq[key] : reverseComplementRefSeq;
}

export function referenceSequenceIncludes(seq, key = '') {
  return forRefSeqIncludes(seq, key) || revRefSeqIncludes(seq, key);
}

const forRefSeqIncludes = memoize((seq, key = '') => {
  return key ? refSeq.ref_seq[key].includes(seq) : refSeq.ref_seq.includes(seq);
});
const revRefSeqIncludes = memoize((seq, key = '') => {
  return key
    ? reverseComplementRefSeq[key].includes(seq)
    : reverseComplementRefSeq.includes(seq);
});

export const queryReferenceSequence = memoize((seq, key = '') => {
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
      res = lookForSeq(newseq, key);
    }
    return res;
  } else {
    // No degenerate nucleotides, keep original logic
    return lookForSeq(seq, key);
  }
});

const forRefSeqIndexOf = memoize((seq, key = '') => {
  return key ? refSeq.ref_seq[key].indexOf(seq) : refSeq.ref_seq.indexOf(seq);
});
const revRefSeqIndexOf = memoize((seq, key = '') => {
  return key
    ? reverseComplementRefSeq[key].indexOf(seq)
    : reverseComplementRefSeq.indexOf(seq);
});

function replaceChar(str, ind, char) {
  return str.substr(0, ind) + char + str.substr(ind + 1);
}

const lookForSeq = memoize((seq, key = '') => {
  if (forRefSeqIncludes(seq, key)) {
    return [
      forRefSeqIndexOf(seq, key) + 1,
      forRefSeqIndexOf(seq, key) + seq.length,
    ];
  } else if (revRefSeqIncludes(seq, key)) {
    return key
      ? [
          refSeq.ref_seq[key].length -
            revRefSeqIndexOf(seq, key) -
            seq.length +
            1,
          refSeq.ref_seq[key].length - revRefSeqIndexOf(seq, key),
        ]
      : [
          refSeq.ref_seq.length - revRefSeqIndexOf(seq) - seq.length + 1,
          refSeq.ref_seq.length - revRefSeqIndexOf(seq),
        ];
  } else {
    return 0;
  }
});

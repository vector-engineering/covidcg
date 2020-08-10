import refSeq from '../../static_data/reference.json';

import _ from 'underscore';
import { reverseComplement } from './string';

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

const forRefSeqIncludes = _.memoize((seq) => {
  return refSeq.ref_seq.includes(seq);
});
const revRefSeqIncludes = _.memoize((seq) => {
  return reverseComplementRefSeq.includes(seq);
});

export const queryReferenceSequence = _.memoize((seq) => {
  if (forRefSeqIncludes(seq)) {
    return [forRefSeqIndexOf(seq) + 1, forRefSeqIndexOf(seq) + seq.length];
  } else if (revRefSeqIncludes(seq)) {
    return [
      refSeq.ref_seq.length - revRefSeqIndexOf(seq) - seq.length + 1,
      refSeq.ref_seq.length - revRefSeqIndexOf(seq),
    ];
  }
});

const forRefSeqIndexOf = _.memoize((seq) => {
  return refSeq.ref_seq.indexOf(seq);
});
const revRefSeqIndexOf = _.memoize((seq) => {
  return reverseComplementRefSeq.indexOf(seq);
});

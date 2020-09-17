/*
 * Load SNP data, and map integers -> SNP strings
 */

// import dnaSnvMap from '../../data/dna_snp_map.json';
// import geneAaSnvMap from '../../data/gene_aa_snp_map.json';
// import proteinAaSnvMap from '../../data/protein_aa_snp_map.json';

const dnaSnvMap = {};
const geneAaSnvMap = {};
const proteinAaSnvMap = {};

import _ from 'underscore';

import { getGene, getProtein } from './gene_protein';

import { snpColorArray } from '../constants/colors';
import { GROUPS } from '../constants/groups';
import { DNA_OR_AA } from '../constants/config';

// Make a SNV -> color map

let snvColorInd = 0;
const _getSnvColor = _.memoize(() => {
  const color = snpColorArray[snvColorInd++];
  if (snvColorInd === snpColorArray.length) {
    snvColorInd = 0;
  }
  return color;
});

const snvColorMap = {};
snvColorMap[GROUPS.REFERENCE_GROUP] = _getSnvColor('Reference');
snvColorMap[GROUPS.OTHER_GROUP] = '#AAA';
snvColorMap[GROUPS.NONE_GROUP] = '#AAA';
snvColorMap[GROUPS.ALL_OTHER_GROUP] = '#AAA';
Object.keys(dnaSnvMap).forEach((snv) => {
  snvColorMap[snv] = _getSnvColor(snv);
});
Object.keys(geneAaSnvMap).forEach((snv) => {
  snvColorMap[snv] = _getSnvColor(snv);
});
Object.keys(proteinAaSnvMap).forEach((snv) => {
  snvColorMap[snv] = _getSnvColor(snv);
});

export const getSnvColor = (snv) => {
  return snvColorMap[snv];
};

// export function loadAaSnpMap() {
//   return aaSnpMap;
// }
// export function loadDnaSnpMap() {
//   return dnaSnpMap;
// }

// Re-map so it's integer -> SNP
let intToDnaSnvMap = {};
let intToGeneAaSnvMap = {};
let intToProteinAaSnvMap = {};

let snvId, split, aaRangeInd;

Object.keys(dnaSnvMap).forEach((snv) => {
  snvId = parseInt(dnaSnvMap[snv]);
  intToDnaSnvMap[snvId] = {};
  // Store the entire SNV string
  intToDnaSnvMap[snvId]['snp_str'] = snv;

  // Each SNV is broken up by pos|ref|alt
  split = snv.split('|');

  // Store all the parts
  // Positions are 1-indexed
  intToDnaSnvMap[snvId]['pos'] = parseInt(split[0]);
  intToDnaSnvMap[snvId]['ref'] = split[1];
  intToDnaSnvMap[snvId]['alt'] = split[2];
});
Object.keys(geneAaSnvMap).forEach((snv) => {
  snvId = parseInt(geneAaSnvMap[snv]);
  intToGeneAaSnvMap[snvId] = {};
  // Store the entire SNV string
  intToGeneAaSnvMap[snvId]['snp_str'] = snv;

  // Each SNV is broken up by gene|pos|ref|alt
  split = snv.split('|');

  // Store all the parts
  intToGeneAaSnvMap[snvId]['gene'] = split[0];
  intToGeneAaSnvMap[snvId]['pos'] = parseInt(split[1]);
  intToGeneAaSnvMap[snvId]['ref'] = split[2];
  intToGeneAaSnvMap[snvId]['alt'] = split[3];

  // Get coordinates in NT (from start of codon)
  aaRangeInd = getGene(split[0]).aa_ranges.reduce((_aaRangeInd, range, ind) => {
    return intToGeneAaSnvMap[snvId]['pos'] >= range[0] &&
      intToGeneAaSnvMap[snvId]['pos'] <= range[1]
      ? ind
      : _aaRangeInd;
  }, 0);
  intToGeneAaSnvMap[snvId]['nt_pos'] =
    getGene(split[0]).ranges[aaRangeInd][0] +
    (intToGeneAaSnvMap[snvId]['pos'] -
      getGene(split[0]).aa_ranges[aaRangeInd][0]) *
      3;
});
Object.keys(proteinAaSnvMap).forEach((snv) => {
  snvId = parseInt(proteinAaSnvMap[snv]);
  intToProteinAaSnvMap[snvId] = {};
  // Store the entire SNV string
  intToProteinAaSnvMap[snvId]['snp_str'] = snv;

  // Each SNV is broken up by gene|pos|ref|alt
  split = snv.split('|');

  // Store all the parts
  intToProteinAaSnvMap[snvId]['protein'] = split[0];
  intToProteinAaSnvMap[snvId]['pos'] = parseInt(split[1]);
  intToProteinAaSnvMap[snvId]['ref'] = split[2];
  intToProteinAaSnvMap[snvId]['alt'] = split[3];

  // Get coordinates in NT (from start of codon)
  aaRangeInd = getProtein(split[0]).aa_ranges.reduce(
    (_aaRangeInd, range, ind) => {
      return intToProteinAaSnvMap[snvId]['pos'] >= range[0] &&
        intToProteinAaSnvMap[snvId]['pos'] <= range[1]
        ? ind
        : _aaRangeInd;
    },
    0
  );
  intToProteinAaSnvMap[snvId]['nt_pos'] =
    getProtein(split[0]).ranges[aaRangeInd][0] +
    (intToProteinAaSnvMap[snvId]['pos'] -
      getProtein(split[0]).aa_ranges[aaRangeInd][0]) *
      3;
});

// console.log(intToDnaSnvMap);
// console.log(intToGeneAaSnvMap);
// console.log(intToProteinAaSnvMap);

export function intToDnaSnv(dnaSnvId) {
  if (dnaSnvId === -1) {
    return { snp_str: GROUPS.REFERENCE_GROUP };
  }
  return intToDnaSnvMap[dnaSnvId];
}
export function intToGeneAaSnv(aaSnvId) {
  if (aaSnvId === -1) {
    return { snp_str: GROUPS.REFERENCE_GROUP };
  }
  return intToGeneAaSnvMap[aaSnvId];
}
export function intToProteinAaSnv(aaSnvId) {
  if (aaSnvId === -1) {
    return { snp_str: GROUPS.REFERENCE_GROUP };
  }
  return intToProteinAaSnvMap[aaSnvId];
}

export function dnaSnvToInt(dnaSnv) {
  return dnaSnvMap[dnaSnv];
}
export function geneAaSnvToInt(geneAaSnv) {
  return geneAaSnvMap[geneAaSnv];
}
export function proteinAaSnvToInt(proteinAaSnv) {
  return proteinAaSnvMap[proteinAaSnv];
}

export function formatSnv(snvStr, dnaOrAa) {
  // Don't do this if it's a special group
  if (Object.values(GROUPS).includes(snvStr)) {
    return snvStr;
  }

  // Print as REF POS ALT
  // i.e., 23403|A|G -> A23403G, S|614|D|G -> S · D614G
  const chunks = snvStr.split('|');
  if (dnaOrAa === DNA_OR_AA.DNA) {
    return `${chunks[1]}${chunks[0]}${chunks[2]}`;
  } else if (dnaOrAa === DNA_OR_AA.AA) {
    return `${chunks[0]} · ${chunks[2]}${chunks[1]}${chunks[3]}`;
  }
}

/*
 * Load SNP data, and map integers -> SNP strings
 */

import dnaSnpMap from '../../data/dna_snp_map.json';
import geneAaSnpMap from '../../data/gene_aa_snp_map.json';
import proteinAaSnpMap from '../../data/protein_aa_snp_map.json';

import _ from 'underscore';

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
Object.keys(dnaSnpMap).forEach((snv) => {
  snvColorMap[snv] = _getSnvColor(snv);
});
Object.keys(geneAaSnpMap).forEach((snv) => {
  snvColorMap[snv] = _getSnvColor(snv);
});
Object.keys(proteinAaSnpMap).forEach((snv) => {
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
let intToDnaSnpMap = {};
let intToGeneAaSnpMap = {};
let intToProteinAaSnpMap = {};

let snpId = -1;
let split = [];

Object.keys(dnaSnpMap).forEach((snp) => {
  snpId = parseInt(dnaSnpMap[snp]);
  intToDnaSnpMap[snpId] = {};
  // Store the entire SNP string
  intToDnaSnpMap[snpId]['snp_str'] = snp;

  // Each SNP is broken up by pos|ref|alt
  split = snp.split('|');

  // Store all the parts
  // Positions are 1-indexed
  intToDnaSnpMap[snpId]['pos'] = parseInt(split[0]);
  intToDnaSnpMap[snpId]['ref'] = split[1];
  intToDnaSnpMap[snpId]['alt'] = split[2];
});
Object.keys(geneAaSnpMap).forEach((snp) => {
  snpId = parseInt(geneAaSnpMap[snp]);
  intToGeneAaSnpMap[snpId] = {};
  // Store the entire SNP string
  intToGeneAaSnpMap[snpId]['snp_str'] = snp;

  // Each SNP is broken up by gene|pos|ref|alt
  split = snp.split('|');

  // Store all the parts
  intToGeneAaSnpMap[snpId]['gene'] = split[0];
  intToGeneAaSnpMap[snpId]['pos'] = parseInt(split[1]);
  intToGeneAaSnpMap[snpId]['ref'] = split[2];
  intToGeneAaSnpMap[snpId]['alt'] = split[3];
});
Object.keys(proteinAaSnpMap).forEach((snp) => {
  snpId = parseInt(proteinAaSnpMap[snp]);
  intToProteinAaSnpMap[snpId] = {};
  // Store the entire SNP string
  intToProteinAaSnpMap[snpId]['snp_str'] = snp;

  // Each SNP is broken up by gene|pos|ref|alt
  split = snp.split('|');

  // Store all the parts
  intToProteinAaSnpMap[snpId]['protein'] = split[0];
  intToProteinAaSnpMap[snpId]['pos'] = parseInt(split[1]);
  intToProteinAaSnpMap[snpId]['ref'] = split[2];
  intToProteinAaSnpMap[snpId]['alt'] = split[3];
});

// console.log(intToDnaSnpMap);
// console.log(intToGeneAaSnpMap);
// console.log(intToProteinAaSnpMap);

export function intToDnaSnp(dnaSnpId) {
  if (dnaSnpId === -1) {
    return { snp_str: GROUPS.REFERENCE_GROUP };
  }
  return intToDnaSnpMap[dnaSnpId];
}
export function intToGeneAaSnp(aaSnpId) {
  if (aaSnpId === -1) {
    return { snp_str: GROUPS.REFERENCE_GROUP };
  }
  return intToGeneAaSnpMap[aaSnpId];
}
export function intToProteinAaSnp(aaSnpId) {
  if (aaSnpId === -1) {
    return { snp_str: GROUPS.REFERENCE_GROUP };
  }
  return intToProteinAaSnpMap[aaSnpId];
}

export function dnaSnpToInt(dnaSnp) {
  return dnaSnpMap[dnaSnp];
}
export function geneAaSnpToInt(geneAaSnp) {
  return geneAaSnpMap[geneAaSnp];
}
export function proteinAaSnpToInt(proteinAaSnp) {
  return proteinAaSnpMap[proteinAaSnp];
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

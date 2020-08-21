/*
 * Load SNP data, and map integers -> SNP strings
 */

import geneAaSnpMap from '../../data/gene_aa_snp_map.json';
import proteinAaSnpMap from '../../data/protein_aa_snp_map.json';
import dnaSnpMap from '../../data/dna_snp_map.json';

import { REFERENCE_GROUP } from '../constants/groups';

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
    return { snp_str: REFERENCE_GROUP };
  }
  return intToDnaSnpMap[dnaSnpId];
}
export function intToGeneAaSnp(aaSnpId) {
  if (aaSnpId === -1) {
    return { snp_str: REFERENCE_GROUP };
  }
  return intToGeneAaSnpMap[aaSnpId];
}
export function intToProteinAaSnp(aaSnpId) {
  if (aaSnpId === -1) {
    return { snp_str: REFERENCE_GROUP };
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

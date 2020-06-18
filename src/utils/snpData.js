/*
 * Load SNP data, and map integers -> SNP strings
 */

import aaSnpMap from '../../data/aa_snp_map.json';
import dnaSnpMap from '../../data/dna_snp_map.json';

import { getAllGenes } from './gene';

// export function loadAaSnpMap() {
//   return aaSnpMap;
// }
// export function loadDnaSnpMap() {
//   return dnaSnpMap;
// }

// Re-map so it's integer -> SNP
let intToDnaSnpMap = {};
let intToAaSnpMap = {};

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
  intToDnaSnpMap[snpId]['pos'] = parseInt(split[0]);
  intToDnaSnpMap[snpId]['ref'] = split[1];
  intToDnaSnpMap[snpId]['alt'] = split[2];
});
Object.keys(aaSnpMap).forEach((snp) => {
  snpId = parseInt(aaSnpMap[snp]);
  intToAaSnpMap[snpId] = {};
  // Store the entire SNP string
  intToAaSnpMap[snpId]['snp_str'] = snp;

  // Each SNP is broken up by gene|pos|ref|alt
  split = snp.split('|');

  // Store all the parts
  intToAaSnpMap[snpId]['gene'] = split[0];
  intToAaSnpMap[snpId]['pos'] = parseInt(split[1]);
  intToAaSnpMap[snpId]['ref'] = split[2];
  intToAaSnpMap[snpId]['alt'] = split[3];
});

// console.log(intToDnaSnpMap);
// console.log(intToAaSnpMap);

export function intToDnaSnp(dnaSnpId) {
  if (dnaSnpId === -1) {
    return { snp_str: 'Reference' };
  }

  return intToDnaSnpMap[dnaSnpId];
}
export function intToAaSnp(aaSnpId) {
  if (aaSnpId === -1) {
    return { snp_str: 'Reference' };
  }

  return intToAaSnpMap[aaSnpId];
}

// Build a map of snp_ids -> gene
let geneToDnaSnp = {};
let geneToAaSnp = {};
let allGenes = getAllGenes();

let snpObj = {};
allGenes.forEach((geneObj) => {
  geneToDnaSnp[geneObj.gene] = [];
  geneToAaSnp[geneObj.gene] = [];

  Object.keys(intToDnaSnpMap).forEach((dnaSnpId) => {
    snpObj = intToDnaSnpMap[dnaSnpId];

    // Check that the SNP position is within the boundaries for this gene
    if (snpObj.pos >= geneObj.start && snpObj.pos <= geneObj.end) {
      geneToDnaSnp[geneObj.gene].push(parseInt(dnaSnpId));
    }
  });

  Object.keys(intToAaSnpMap).forEach((aaSnpId) => {
    snpObj = intToAaSnpMap[aaSnpId];

    // Check that the AA SNP gene is the same as this gene
    if (snpObj.gene === geneObj.gene) {
      geneToAaSnp[geneObj.gene].push(parseInt(aaSnpId));
    }
  });
});

// console.log(geneToDnaSnp);
// console.log(geneToAaSnp);

// Is this SNP ID within this gene?
export function dnaSnpInGene(dnaSnpId, geneName) {
  return geneToDnaSnp[geneName].includes(dnaSnpId);
}
export function aaSnpInGene(aaSnpId, geneName) {
  return intToAaSnpMap[aaSnpId]['gene'] === geneName;
}

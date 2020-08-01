/*
 * Load global counts per group (lineage, clade, SNV, etc)
 */

import globalGroupCounts from '../../data/global_group_counts.json';
import { intToDnaSnp, intToGeneAaSnp, intToProteinAaSnp } from './snpData';

// Make a copy
const processedGlobalGroupCounts = Object.assign({}, globalGroupCounts);

// Replace integer IDs with SNP strings
let snpString;
Object.keys(processedGlobalGroupCounts.dna_snp).forEach((snpId) => {
  snpString = intToDnaSnp(snpId).snp_str;
  processedGlobalGroupCounts.dna_snp[snpString] =
    processedGlobalGroupCounts.dna_snp[snpId];
  delete processedGlobalGroupCounts.dna_snp[snpId];
});

Object.keys(processedGlobalGroupCounts.gene_aa_snp).forEach((snpId) => {
  snpString = intToGeneAaSnp(snpId).snp_str;
  processedGlobalGroupCounts.gene_aa_snp[snpString] =
    processedGlobalGroupCounts.gene_aa_snp[snpId];
  delete processedGlobalGroupCounts.gene_aa_snp[snpId];
});

Object.keys(processedGlobalGroupCounts.protein_aa_snp).forEach((snpId) => {
  snpString = intToProteinAaSnp(snpId).snp_str;
  processedGlobalGroupCounts.protein_aa_snp[snpString] =
    processedGlobalGroupCounts.protein_aa_snp[snpId];
  delete processedGlobalGroupCounts.protein_aa_snp[snpId];
});

// console.log(processedGlobalGroupCounts);

export function getGlobalGroupCounts() {
  return processedGlobalGroupCounts;
}

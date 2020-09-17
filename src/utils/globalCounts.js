/*
 * Load global counts per group (lineage, clade, SNV, etc)
 */

//import globalGroupCounts from '../../data/global_group_counts.json';
const globalGroupCounts = [];

// Make a copy
const processedGlobalGroupCounts = Object.assign({}, globalGroupCounts);

// Replace integer IDs with SNP strings
Object.keys(processedGlobalGroupCounts.dna_snp).forEach((snpId) => {
  processedGlobalGroupCounts.dna_snp[snpId.toString()] =
    processedGlobalGroupCounts.dna_snp[snpId];
});

Object.keys(processedGlobalGroupCounts.gene_aa_snp).forEach((snpId) => {
  processedGlobalGroupCounts.gene_aa_snp[snpId.toString()] =
    processedGlobalGroupCounts.gene_aa_snp[snpId];
});

Object.keys(processedGlobalGroupCounts.protein_aa_snp).forEach((snpId) => {
  processedGlobalGroupCounts.protein_aa_snp[snpId.toString()] =
    processedGlobalGroupCounts.protein_aa_snp[snpId];
});

// console.log(processedGlobalGroupCounts);

export function getGlobalGroupCounts() {
  return processedGlobalGroupCounts;
}

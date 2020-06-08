/*
 * Load SNP data, and map integers -> SNP strings
 */

import aa_snp_map from '../../data/aa_snp_map.json';
import dna_snp_map from '../../data/dna_snp_map.json';

import { getAllGenes } from './gene';

// export function load_aa_snp_map() {
//   return aa_snp_map;
// }
// export function load_dna_snp_map() {
//   return dna_snp_map;
// }

// Re-map so it's integer -> SNP
let int_to_dna_snp_map = {};
let int_to_aa_snp_map = {};

let snp_id = -1;
let split = [];

Object.keys(dna_snp_map).forEach((snp) => {
  snp_id = parseInt(dna_snp_map[snp]);
  int_to_dna_snp_map[snp_id] = {};
  // Store the entire SNP string
  int_to_dna_snp_map[snp_id]['snp_str'] = snp;

  // Each SNP is broken up by pos|ref|alt
  split = snp.split('|');

  // Store all the parts
  int_to_dna_snp_map[snp_id]['pos'] = parseInt(split[0]);
  int_to_dna_snp_map[snp_id]['ref'] = split[1];
  int_to_dna_snp_map[snp_id]['alt'] = split[2];
});
Object.keys(aa_snp_map).forEach((snp) => {
  snp_id = parseInt(aa_snp_map[snp]);
  int_to_aa_snp_map[snp_id] = {};
  // Store the entire SNP string
  int_to_aa_snp_map[snp_id]['snp_str'] = snp;

  // Each SNP is broken up by gene|pos|ref|alt
  split = snp.split('|');

  // Store all the parts
  int_to_aa_snp_map[snp_id]['gene'] = split[0];
  int_to_aa_snp_map[snp_id]['pos'] = parseInt(split[1]);
  int_to_aa_snp_map[snp_id]['ref'] = split[2];
  int_to_aa_snp_map[snp_id]['alt'] = split[3];
});

// console.log(int_to_dna_snp_map);
// console.log(int_to_aa_snp_map);

export function int_to_dna_snp(dna_snp_id) {
  if (dna_snp_id === -1) {
    return { snp_str: 'Reference' };
  }

  return int_to_dna_snp_map[dna_snp_id];
}
export function int_to_aa_snp(aa_snp_id) {
  if (aa_snp_id === -1) {
    return { snp_str: 'Reference' };
  }

  return int_to_aa_snp_map[aa_snp_id];
}

// Build a map of snp_ids -> gene
let gene_to_dna_snp = {};
let gene_to_aa_snp = {};
let allGenes = getAllGenes();

let snp_obj = {};
allGenes.forEach((gene_obj) => {
  gene_to_dna_snp[gene_obj.gene] = [];
  gene_to_aa_snp[gene_obj.gene] = [];

  Object.keys(int_to_dna_snp_map).forEach((dna_snp_id) => {
    snp_obj = int_to_dna_snp_map[dna_snp_id];

    // Check that the SNP position is within the boundaries for this gene
    if (snp_obj.pos >= gene_obj.start && snp_obj.pos <= gene_obj.end) {
      gene_to_dna_snp[gene_obj.gene].push(parseInt(dna_snp_id));
    }
  });

  Object.keys(int_to_aa_snp_map).forEach((aa_snp_id) => {
    snp_obj = int_to_aa_snp_map[aa_snp_id];

    // Check that the AA SNP gene is the same as this gene
    if (snp_obj.gene === gene_obj.gene) {
      gene_to_aa_snp[gene_obj.gene].push(parseInt(aa_snp_id));
    }
  });
});

// console.log(gene_to_dna_snp);
// console.log(gene_to_aa_snp);

// Is this SNP ID within this gene?
export function dna_snp_in_gene(dna_snp_id, gene_name) {
  return gene_to_dna_snp[gene_name].includes(dna_snp_id);
}
export function aa_snp_in_gene(aa_snp_id, gene_name) {
  return int_to_aa_snp_map[aa_snp_id]['gene'] === gene_name;
}

import lineageSnpData from '../../data/lineage_snp.json';
import cladeSnpData from '../../data/clade_snp.json';
import refSeq from '../../static_data/reference.json';

import _ from 'underscore';
import { intToDnaSnp, intToGeneAaSnp, intToProteinAaSnp } from './snpData';

export function loadLineageSnp() {
  return lineageSnpData;
}
export function loadCladeSnp() {
  return cladeSnpData;
}

export function getReferenceSequence() {
  return refSeq['ref_seq'];
}

function getGroup(groupKey, group) {
  let snpData;
  if (groupKey === 'lineage') {
    snpData = lineageSnpData;
  } else if (groupKey === 'clade') {
    snpData = cladeSnpData;
  }
  let findObj = {};
  findObj[groupKey] = group;

  return _.findWhere(snpData, findObj);
}

export function getDnaSnpsFromGroup(groupKey, group) {
  let groupObj = getGroup(groupKey, group);
  if (groupObj === undefined) {
    return [];
  }

  let snpIds = groupObj.dna_snp_ids;
  snpIds = _.reject(snpIds, (snpId) => snpId === '');
  return _.map(snpIds, (snpId) => intToDnaSnp(parseInt(snpId)));
}

export function getGeneAaSnpsFromGroup(groupKey, group) {
  let groupObj = getGroup(groupKey, group);
  if (groupObj === undefined) {
    return [];
  }
  let snpIds = groupObj.gene_aa_snp_ids;
  snpIds = _.reject(snpIds, (snpId) => snpId === '');
  return _.map(snpIds, (snpId) => intToGeneAaSnp(parseInt(snpId)));
}

export function getProteinAaSnpsFromGroup(groupKey, group) {
  let groupObj = getGroup(groupKey, group);
  if (groupObj === undefined) {
    return [];
  }
  let snpIds = groupObj.protein_aa_snp_ids;
  snpIds = _.reject(snpIds, (snpId) => snpId === '');
  return _.map(snpIds, (snpId) => intToProteinAaSnp(parseInt(snpId)));
}

// export function getLineagesWithDnaSnpInGene(gene) {
//   let validLineages = new Set();

//   // Get lineages whose positions fall within start -- end
//   lineageSnpData.forEach((lineageObj) => {
//     // Get SNP objects from this list of IDs
//     let snps = getDnaSnpsFromLineage(lineageObj.lineage);

//     snps.forEach((snp) => {
//       // If any one of it's positions is in the range, then add this lineage
//       if (snp.pos >= gene.start && snp.pos <= gene.end) {
//         validLineages.add(lineageObj.lineage);
//       }
//     });
//   });

//   return Array.from(validLineages);
// }

// export function getLineagesWithAaSnpInGene(gene) {
//   let lineageData = loadLineageSnp();
//   let validLineages = new Set();

//   // Get lineages whose positions fall within start -- end
//   lineageData.forEach((lineageObj) => {
//     // Get SNP objects from this list of IDs
//     let snps = getGeneAaSnpsFromLineage(lineageObj.lineage);

//     snps.forEach((snp) => {
//       // If the SNP is on the gene of interest, then add it
//       if (snp.gene === gene.gene) {
//         validLineages.add(lineageObj.lineage);
//       }
//     });
//   });

//   return Array.from(validLineages);
// }

// export function getLineagesWithAaSnpInProtein(gene) {
//   let lineageData = loadLineageSnp();
//   let validLineages = new Set();

//   // Get lineages whose positions fall within start -- end
//   lineageData.forEach((lineageObj) => {
//     // Get SNP objects from this list of IDs
//     let snps = getProteinAaSnpsFromLineage(lineageObj.lineage);

//     snps.forEach((snp) => {
//       // If the SNP is on the gene of interest, then add it
//       if (snp.gene === gene.gene) {
//         validLineages.add(lineageObj.lineage);
//       }
//     });
//   });

//   return Array.from(validLineages);
// }

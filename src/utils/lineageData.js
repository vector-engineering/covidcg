import lineageSnpData from '../../data/lineage_snp.json';
import refSeq from '../../static_data/reference.json';

import _ from 'underscore';
import { intToDnaSnp, intToAaSnp } from './snpData';

export function loadLineageSnp() {
  return lineageSnpData;
}

export function getReferenceSequence() {
  return refSeq['ref_seq'];
}

export function getDnaSnpsFromLineage(lineage) {
  let lineageObj = _.findWhere(lineageSnpData, { lineage: lineage });
  if (lineageObj === undefined) {
    return [];
  }
  let snpIds = lineageObj.dna_snp_ids;
  snpIds = _.reject(snpIds, (snpId) => snpId === '');
  return _.map(snpIds, (snpId) => intToDnaSnp(parseInt(snpId)));
}

export function getAaSnpsFromLineage(lineage) {
  let lineageObj = _.findWhere(lineageSnpData, { lineage: lineage });
  if (lineageObj === undefined) {
    return [];
  }
  let snpIds = lineageObj.aa_snp_ids;
  snpIds = _.reject(snpIds, (snpId) => snpId === '');
  return _.map(snpIds, (snpId) => intToAaSnp(parseInt(snpId)));
}

export function getLineagesWithDnaSnpInGene(gene) {
  let validLineages = new Set();

  // Get lineages whose positions fall within start -- end
  lineageSnpData.forEach((lineageObj) => {
    // Get SNP objects from this list of IDs
    let snps = getDnaSnpsFromLineage(lineageObj.lineage);

    snps.forEach((snp) => {
      // If any one of it's positions is in the range, then add this lineage
      if (snp.pos >= gene.start && snp.pos <= gene.end) {
        validLineages.add(lineageObj.lineage);
      }
    });
  });

  return Array.from(validLineages);
}

export function getLineagesWithAaSnpInGene(gene) {
  let lineageData = loadLineageSnp();
  let validLineages = new Set();

  // Get lineages whose positions fall within start -- end
  lineageData.forEach((lineageObj) => {
    // Get SNP objects from this list of IDs
    let snps = getAaSnpsFromLineage(lineageObj.lineage);

    snps.forEach((snp) => {
      // If the SNP is on the gene of interest, then add it
      if (snp.gene === gene.gene) {
        validLineages.add(lineageObj.lineage);
      }
    });
  });

  return Array.from(validLineages);
}

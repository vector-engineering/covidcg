import lineageData from '../../data/lineage_dna_snp.json';
import refSeq from '../../data/reference.json';
import _ from 'underscore';

export function loadLineageData() {
  return lineageData;
}

export function getReferenceSequence() {
  return refSeq['ref_seq'];
}

export function getLineagesFromGene(gene) {
  let start_pos = gene.start;
  let end_pos = gene.end;

  let lineageData = loadLineageData();
  let validLineages = new Set();

  // Get lineages whose positions fall within start -- end
  lineageData.forEach(row => {
    // If any one of it's positions is in the range, then add it
    if(row.pos >= start_pos && row.pos <= end_pos) {
      validLineages.add(row.lineage);
    }
  });
  
  return Array.from(validLineages);
}
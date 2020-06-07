import lineageDnaSnp from '../../data/lineage_dna_snp.json';
import lineageAaSnp from '../../data/lineage_aa_snp.json';
import refSeq from '../../data/reference.json';

export function loadLineageDnaSnp() {
  return lineageDnaSnp;
}
export function loadLineageAaSnp() {
  return lineageAaSnp;
}

export function getReferenceSequence() {
  return refSeq['ref_seq'];
}

export function getLineagesWithDnaSnpInGene(gene) {
  let lineageData = loadLineageDnaSnp();
  let validLineages = new Set();

  // Get lineages whose positions fall within start -- end
  lineageData.forEach((row) => {
    // If any one of it's positions is in the range, then add it
    if (row.pos >= gene.start && row.pos <= gene.end) {
      validLineages.add(row.lineage);
    }
  });

  return Array.from(validLineages);
}

export function getLineagesWithAaSnpInGene(gene) {
  let lineageData = loadLineageAaSnp();
  let validLineages = new Set();

  // Get lineages whose positions fall within start -- end
  lineageData.forEach((row) => {
    // If the SNP is on the gene of interest, then add it
    if (row.gene === gene.gene) {
      validLineages.add(row.lineage);
    }
  });

  return Array.from(validLineages);
}

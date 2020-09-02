import genes from '../../static_data/genes.json';
import proteins from '../../static_data/proteins.json';

let processedGenes = genes;

// Parse segments
processedGenes = processedGenes.map((gene) => {
  gene.ranges = gene.segments.split(';').map((segment) => {
    return segment.split('..').map((pos) => parseInt(pos));
  });
  let curResidueIndex = 1;
  gene.aa_ranges = gene.protein_coding
    ? gene.ranges.map((range) => {
        const aa_range = [curResidueIndex, (range[1] - range[0] + 1) / 3];
        curResidueIndex = aa_range[1] + 1;
        return aa_range;
      })
    : null;
  gene.len_nt = gene.ranges.reduce((len, range) => {
    return len + (range[1] - range[0]) + 1;
  }, 0);
  gene.len_aa = gene.protein_coding ? Math.floor(gene.len_nt / 3) : null;
  return gene;
});

export function getAllGenes() {
  return processedGenes;
}

const geneMap = {
  'All Genes': {
    gene: 'All Genes',
    ranges: [[1, 30000]],
    domains: [],
  },
};
processedGenes.forEach((gene) => {
  geneMap[gene.gene] = gene;
});

export function getGene(gene) {
  // Get the selected gene object
  return geneMap[gene];
}

let processedProteins = proteins;

// Parse segments
processedProteins = processedProteins.map((protein) => {
  protein.ranges = protein.segments.split(';').map((segment) => {
    return segment.split('..').map((pos) => parseInt(pos));
  });
  let curResidueIndex = 1;
  protein.aa_ranges = protein.ranges.map((range) => {
    const aa_range = [curResidueIndex, (range[1] - range[0] + 1) / 3];
    curResidueIndex = aa_range[1] + 1;
    return aa_range;
  });
  protein.len_nt = protein.ranges.reduce((len, range) => {
    return len + (range[1] - range[0]) + 1;
  }, 0);
  protein.len_aa = Math.floor(protein.len_nt / 3);
  return protein;
});

// console.log(processedProteins);

export function getAllProteins() {
  return processedProteins;
}

const proteinMap = {
  'All Proteins': {
    protein: 'All Proteins',
    ranges: [[1, 30000]],
    domains: [],
  },
};
processedProteins.forEach((protein) => {
  proteinMap[protein.protein] = protein;
});

export function getProtein(_protein) {
  return proteinMap[_protein];
}

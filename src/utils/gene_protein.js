//__VIRUS__ is replaced with config.virus at build time
import genes from '../../static_data/__VIRUS__/genes_processed.json';
import proteins from '../../static_data/__VIRUS__/proteins_processed.json';

/* function processFeatures(features) {
  return features.map((feature) => {
    feature.ranges = feature.segments;
    let curResidueIndex = 1;
    feature.aa_ranges = feature.protein_coding
      ? feature.ranges.map((range) => {
          const aa_range = [
            curResidueIndex,
            curResidueIndex - 1 + (range[1] - range[0] + 1) / 3,
          ];
          curResidueIndex = aa_range[1] + 1;
          return aa_range;
        })
      : null;
    feature.len_nt = feature.ranges.reduce((len, range) => {
      return len + (range[1] - range[0]) + 1;
    }, 0);
    feature.len_aa = feature.protein_coding
      ? Math.floor(feature.len_nt / 3)
      : null;
    return feature;
  });
} */

export function getAllGenes() {
  return genes;
}

export const geneMap = {
  'All Genes': {
    name: 'All Genes',
    ranges: [[1, 30000]],
    domains: [],
  },
};
genes.forEach((gene) => {
  geneMap[gene.name] = gene;
});

export function getGene(gene) {
  // Get the selected gene object
  return geneMap[gene];
}

export function getAllProteins() {
  return proteins;
}

export const proteinMap = {
  'All Proteins': {
    name: 'All Proteins',
    ranges: [[1, 30000]],
    domains: [],
  },
};
proteins.forEach((protein) => {
  proteinMap[protein.name] = protein;
});

export function getProtein(_protein) {
  return proteinMap[_protein];
}

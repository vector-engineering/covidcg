<<<<<<< HEAD
import genes from '../../static_data/genes.json';
import proteins from '../../static_data/proteins.json';

function processFeatures(features) {
=======
import genes from '../../static_data/genes_processed.json';
import proteins from '../../static_data/proteins_processed.json';

/* function processFeatures(features) {
>>>>>>> Updated genes and proteins json to work with new entropyplot
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

<<<<<<< HEAD
let processedGenes = processFeatures(genes);
=======
let processedGenes = genes;
>>>>>>> Updated genes and proteins json to work with new entropyplot

export function getAllGenes() {
  return processedGenes;
}

export const geneMap = {
  'All Genes': {
    name: 'All Genes',
    ranges: [[1, 30000]],
    domains: [],
  },
};
processedGenes.forEach((gene) => {
  geneMap[gene.name] = gene;
});

export function getGene(gene) {
  // Get the selected gene object
  return geneMap[gene];
}

<<<<<<< HEAD
let processedProteins = processFeatures(proteins);
=======
let processedProteins = proteins;
>>>>>>> Updated genes and proteins json to work with new entropyplot

export function getAllProteins() {
  return processedProteins;
}

export const proteinMap = {
  'All Proteins': {
    name: 'All Proteins',
    ranges: [[1, 30000]],
    domains: [],
  },
};
processedProteins.forEach((protein) => {
  proteinMap[protein.name] = protein;
});

export function getProtein(_protein) {
  return proteinMap[_protein];
}

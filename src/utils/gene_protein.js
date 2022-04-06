//__VIRUS__ is replaced with config.virus at build time
import genes from '../../static_data/__VIRUS__/genes_processed.json';
import proteins from '../../static_data/__VIRUS__/proteins_processed.json';

import { config } from '../config';

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

export function getAllGenes(ref = 'A') {
  if (config.virus === 'sars2') {
    return genes;
  } else {
    return ref === 'A' ? genes['A'] : genes['B'];
  }
}

export const geneMap = {
  'All Genes': {
    name: 'All Genes',
    ranges: [[1, 30000]],
    domains: [],
  },
};

export function getGene(gene, ref = 'A') {
  // Get the selected gene object
  if (!(gene in geneMap)) {
    if (config.virus === 'sars2') {
      genes.forEach((gene) => {
        geneMap[gene.name] = gene;
      });
    } else {
      if (!(ref in geneMap)) {
        const selectedGenes = ref === 'A' ? genes['A'] : genes['B'];
        geneMap[ref] = {};
        selectedGenes.forEach((gene) => {
          geneMap[ref][gene.name] = gene;
        });
      }
    }
  }
  return config.virus === 'sars2' ? geneMap[gene] : geneMap[ref][gene];
}

export function getAllProteins(ref = 'A') {
  if (config.virus === 'sars2') {
    return proteins;
  } else {
    return ref === 'A' ? proteins['A'] : proteins['B'];
  }
}

export const proteinMap = {
  'All Proteins': {
    name: 'All Proteins',
    ranges: [[1, 30000]],
    domains: [],
  },
};

export function getProtein(_protein, ref = 'A') {
  // Get the selected gene object
  if (!(_protein in proteinMap)) {
    if (config.virus === 'sars2') {
      proteins.forEach((protein) => {
        proteinMap[protein.name] = protein;
      });
    } else {
      if (!(ref in proteinMap)) {
        if (config.virus === 'sars2') {
          proteins.forEach((protein) => {
            proteinMap[protein.name] = protein;
          });
        } else {
          const selectedProteins = ref === 'A' ? proteins['A'] : proteins['B'];
          proteinMap[ref] = {};
          selectedProteins.forEach((protein) => {
            proteinMap[ref][protein.name] = protein;
          });
        }
      }
    }
  }
  return config.virus === 'sars2'
    ? proteinMap[_protein]
    : proteinMap[ref][_protein];
}

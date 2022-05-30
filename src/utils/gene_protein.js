//__VIRUS__ is replaced with config.virus at build time
// eslint-disable-next-line import/no-unresolved
import genes from '../../static_data/__VIRUS__/genes_processed.json';
// eslint-disable-next-line import/no-unresolved
import proteins from '../../static_data/__VIRUS__/proteins_processed.json';

// import { config } from '../config';

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

export function getAllGenes(ref) {
  return genes[ref];
}

export function getGene(gene, ref) {
  // Get the selected gene object
  return genes[ref].find((g) => g.name === gene);
}

export function getAllProteins(ref) {
  return proteins[ref];
}

export function getProtein(protein, ref) {
  return proteins[ref].find((p) => p.name === protein);
}

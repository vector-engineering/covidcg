//__VIRUS__ is replaced with config.virus at build time
// eslint-disable-next-line import/no-unresolved
import genes from '../../static_data/__VIRUS__/genes_processed.json';
// eslint-disable-next-line import/no-unresolved
import proteins from '../../static_data/__VIRUS__/proteins_processed.json';

// import { config } from '../config';
import { configStore as initialValues } from '../constants/initialValues';

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

function _getGene(gene, ref) {
  return genes[ref].find((g) => g.name === gene);
}
export function getGene(gene, ref) {
  // Get the selected gene object
  let geneObj = _getGene(gene, ref);
  // If the gene doesn't exist for this reference, try the initial gene
  if (geneObj === undefined) {
    geneObj = _getGene(initialValues.selectedGene.name, ref);
  }
  // If that doesn't exist, then just use the first gene in the list
  // for this reference
  if (geneObj === undefined) {
    geneObj = getAllGenes(ref)[0];
  }

  return geneObj;
}

export function getAllProteins(ref) {
  return proteins[ref];
}

function _getProtein(protein, ref) {
  return proteins[ref].find((p) => p.name === protein);
}
export function getProtein(protein, ref) {
  // Get the selected protein object
  let proteinObj = _getProtein(protein, ref);
  // If the protein doesn't exist for this reference, try the initial protein
  if (proteinObj === undefined) {
    proteinObj = _getProtein(initialValues.selectedProtein.name, ref);
  }
  // If that doesn't exist, then just use the first protein in the list
  // for this reference
  if (proteinObj === undefined) {
    proteinObj = getAllProteins(ref)[0];
  }

  return proteinObj;
}

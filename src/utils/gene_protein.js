//__VIRUS__ is replaced with config.virus at build time
// eslint-disable-next-line import/no-unresolved
import genes from '../../static_data/__VIRUS__/genes_processed.json';
// eslint-disable-next-line import/no-unresolved
import proteins from '../../static_data/__VIRUS__/proteins_processed.json';

// import { config } from '../config';
import { configStore as initialValues } from '../constants/initialValues';

export function getAllGenes(reference) {
  return genes[reference];
}

function _getGene(gene, reference) {
  return genes[reference].find((g) => g.name === gene);
}
export function getGene(gene, reference) {
  // Get the selected gene object
  let geneObj = _getGene(gene, reference);
  // If the gene doesn't exist for this reference, try the initial gene
  if (geneObj === undefined) {
    geneObj = _getGene(initialValues.selectedGene.name, reference);
  }
  // If that doesn't exist, then just use the first gene in the list
  // for this reference
  if (geneObj === undefined) {
    geneObj = getAllGenes(reference)[0];
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

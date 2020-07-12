import _ from 'underscore';
import genes from '../../static_data/genes.json';

export function getAllGenes() {
  return genes;
}

export function getGene(value) {
  // Get the selected gene object
  return (
    // There's no object for "All Genes", so make one now
    value === 'all'
      ? {
          gene: 'All Genes',
          start: 1,
          // Find the largest end position out of all genes
          end: _.reduce(genes, (memo, opt) => Math.max(memo, opt.end), 1),
        }
      : _.findWhere(genes, { gene: value })
  );
}

export function loadGeneOptions() {
  // Load genes and create option for each gene
  let options = [];

  // All Genes option
  options.push({
    label: 'All Genes',
    value: 'all',
  });

  genes.forEach((gene) => {
    options.push({
      label: gene.gene,
      value: gene.gene,
    });
  });

  return options;
}

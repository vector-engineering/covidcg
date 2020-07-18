import _ from 'underscore';
import genes from '../../static_data/genes.json';

export function getAllGenes() {
  return genes;
}

export function getGene(value) {
  // Get the selected gene object
  return (
    // There's no object for "All Genes", so make one now
    value === 'All Genes'
      ? {
          gene: 'All Genes',
          start: 1,
          // Find the largest end position out of all genes
          end: _.reduce(genes, (memo, opt) => Math.max(memo, opt.end), 1),
        }
      : _.findWhere(genes, { gene: value })
  );
}

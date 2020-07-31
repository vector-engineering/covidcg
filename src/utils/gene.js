import _ from 'underscore';
import genes from '../../static_data/genes.json';

let processedGenes = genes;

// Parse segments
processedGenes = _.map(processedGenes, (gene) => {
  gene.ranges = _.map(gene.segments.split(';'), (segment) => {
    return _.map(segment.split('..'), (pos) => parseInt(pos));
  });
  return gene;
});

export function getAllGenes() {
  return processedGenes;
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

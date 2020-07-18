import _ from 'underscore';
import proteins from '../../static_data/proteins.json';

let processedProteins = proteins;

// Parse segments
processedProteins = _.map(processedProteins, (protein) => {
  protein.ranges = _.map(protein.segments.split(';'), (segment) => {
    return _.map(segment.split('..'), (pos) => parseInt(pos));
  });
  return protein;
});

// console.log(processedProteins);

export function getAllProteins() {
  return processedProteins;
}

export function getProtein(_protein) {
  // There's no object for "All Genes", so make one now
  return _protein === 'All Proteins'
    ? {
        protein: 'All Proteins',
        ranges: [[1, 30000]],
      }
    : _.findWhere(processedProteins, { protein: _protein });
}

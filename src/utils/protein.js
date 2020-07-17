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
  return _.findWhere(processedProteins, { protein: _protein });
}

import _ from 'underscore';
import primers from '../../static_data/primers.json';

let processedPrimers = primers;
// console.log(processedPrimers);

// Build primer tree - with Institutions as the roots
let primerSelectTree = [];
let instObj = null;
processedPrimers.forEach((primer) => {
  instObj = _.findWhere(primerSelectTree, { value: primer.Institution });
  if (instObj === undefined) {
    instObj = {
      value: primer.Institution,
      label: primer.Institution,
      level: 'group',
      children: [],
    };
    primerSelectTree.push(instObj);
  }
  instObj.children.push({
    value: primer.Name,
    label: primer.Name,
    level: 'individual',
  });
});

// console.log(primerSelectTree);

export function getAllPrimers() {
  return processedPrimers;
}

export function getPrimerSelectTree() {
  return primerSelectTree;
}

export function getPrimerByName(primerName) {
  return _.findWhere(processedPrimers, { Name: primerName });
}

export function getPrimersByGroup(groupName) {
  return _.where(processedPrimers, { Institution: groupName });
}

export function queryPrimers(queryObj) {
  return _.where(processedPrimers, queryObj);
}

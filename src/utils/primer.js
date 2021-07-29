import primers from '../../static_data/primers.json';

let processedPrimers = primers;
// console.log(processedPrimers);

// Build primer tree - with Institutions as the roots
let primerSelectTree = [];
let instObj = null;
processedPrimers.forEach((primer) => {
  instObj = primerSelectTree.find((node) => node.value === primer.Institution);
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
  return processedPrimers.find((primer) => primer.Name === primerName);
}

export function getPrimersByGroup(groupName) {
  return processedPrimers.filter((primer) => primer.Institution === groupName);
}

// Example usage: queryPrimers({ Institution: "CDC", Name: "N1" })
export function queryPrimers(queryObj) {
  return processedPrimers.find((primer) => {
    return Object.keys(queryObj).every((key) => {
      return primer[key] === queryObj[key];
    });
  });
}

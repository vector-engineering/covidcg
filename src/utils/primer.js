import { config } from '../config';

let processedPrimers = null;
let primerSelectTree = [];
let instObj = null;

if (config.virus === 'sars2') {
  import('../../static_data/primers.json').then((mod) => {
    const primers = mod.default;
    if (primers && primers.length > 0) {
      // Build primer tree - with Institutions as the roots
      processedPrimers = primers;
      primers.forEach((primer) => {
        instObj = primerSelectTree.find(
          (node) => node.value === primer.Institution
        );
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
    }
  });
}

// console.log(primerSelectTree);

export function getAllPrimers() {
  return processedPrimers;
}

export function getPrimerSelectTree() {
  return primerSelectTree;
}

export function getPrimerByName(primerName) {
  return processedPrimers
    ? processedPrimers.find((primer) => primer.Name === primerName)
    : null;
}

export function getPrimersByGroup(groupName) {
  return processedPrimers
    ? processedPrimers.filter((primer) => primer.Institution === groupName)
    : null;
}

// Example usage: queryPrimers({ Institution: "CDC", Name: "N1" })
export function queryPrimers(queryObj) {
  return processedPrimers
    ? processedPrimers.find((primer) => {
        return Object.keys(queryObj).every((key) => {
          return primer[key] === queryObj[key];
        });
      })
    : null;
}

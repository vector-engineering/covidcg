import { action } from 'mobx';

import { rootStoreInstance } from './rootStore';

import { updateURLFromParams } from '../utils/updateQueryParam';

import {
  geneMap,
  proteinMap,
  getGene,
  getProtein,
} from '../utils/gene_protein';

export class URLMonitor {
  urlParams = new URLSearchParams(window.location.search);
  // paramsToDisplay is a list of variables that can be set from the url by store
  paramsToDisplay = {
    configStore: ['groupkey', 'dnaOrAa'],
    groupDataStore: ['selectedGroups'],
  };

  init() {
    // Iterate through params for each store
    Object.keys(this.paramsToDisplay).forEach((store) => {
      this.paramsToDisplay[store].forEach((param) => {
        if (this.urlParams.has(param)) {
          const newValArr = this.urlParams.getAll(param);
          // Make sure it's a valid entry
          if (param === 'selectedGene') {
            // newValArr should have the form ['str']
            // Ignore any other values in newValArr
            if (newValArr[0] in geneMap) {
              // If the specified gene is in the geneMap get the gene
              const gene = getGene(newValArr[0]);
              // Update configStore
            }
            this.urlParams.set(param, newValArr[0]);
          } else if (param === 'selectedProtein') {
            // newValArr should have the form ['str']
            // Ignore any other values in newValArr
            if (newValArr[0] in proteinMap) {
              // If the specified protein is in the proteinMap get the protein
              const protein = getProtein(newValArr[0]);
              // Update configStore
            }
            this.urlParams.set(param, newValArr[0]);
          }
        } else {
          // A parameter was supplied in the url that cannot be set from the url
          // Remove that parameter from the url and do nothing
          this.urlParams.delete(param);
        }
      });
    });

    updateURLFromParams(this.urlParams);
  }

  @action
  updateURL(pending) {
    const store = pending[store];
    this.resetURLParams(store);

    Object.keys(pending).forEach((param) => {
      if (this.paramsToDisplay[store].includes(param)) {
        this.urlParams.set(pending[param]);
      }
    });
  }

  resetURLParams(store) {
    this.paramsToDisplay[store].forEach((param) => {
      this.urlParams.delete(param);
    });
  }
}

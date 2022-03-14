import { action } from 'mobx';

import { rootStoreInstance } from './rootStore';

import { config } from '../config';

import { updateURLFromParams } from '../utils/updateQueryParam';
import { queryPrimers } from '../utils/primer';
import { getLocationByNameAndLevel } from '../utils/location';

import {
  geneMap,
  proteinMap,
  getGene,
  getProtein,
} from '../utils/gene_protein';

import { GEO_LEVELS, TABS } from '../constants/defs.json';

export class URLMonitor {
  urlParams = new URLSearchParams(window.location.search);
  // paramsToDisplay is a list of variables that can be set to/from the url
  paramsToDisplay = {
    configStore: [
      'groupkey',
      'dnaOrAa',
      'selectedGene',
      'selectedProtein',
      'selectedPrimers',
      'selectedGroups',
      'customCoordinates',
      'customSequences',
      'residueCoordinates',
      'coordinateMode',
      'startDate',
      'endDate',
      'submStartDate',
      'submEndDate',
    ],
    groupDataStore: ['selectedGroups'],
  };
  pendingChanges = {
    configStore: {},
    groupDataStore: {},
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
              this.pendingChanges.configStore[param] = gene;
            }
            this.urlParams.set(param, newValArr[0]);
          } else if (param === 'selectedProtein') {
            // newValArr should have the form ['str']
            // Ignore any other values in newValArr
            if (newValArr[0] in proteinMap) {
              // If the specified protein is in the proteinMap get the protein
              const protein = getProtein(newValArr[0]);
              this.pendingChanges.configStore[param] = protein;
            }
            this.urlParams.set(param, newValArr[0]);
          } else if (
            param === 'customCoordinates' ||
            param === 'residueCoordinates'
          ) {
            // If coordinates are specified, save them as an array of numbers
            // Coordinates can stay as a string in the URL
            let arr = [];
            let value = newValArr[0].split(',');
            value.forEach((item, i) => {
              value[i] = parseInt(item);

              if (i % 2 === 1) {
                arr.push([value[i - 1], value[i]]);
              }
            });

            this.pendingChanges.configStore[param] = arr;
          } else if (
            param === 'customSequences' ||
            param === 'selectedGroups'
          ) {
            // Store customSequences as an array of strings
            const value = newValArr[0].split(',');
            this.pendingChanges.configStore[param] = value;
          } else if (param === 'selectedPrimers') {
            const value = newValArr[0].split(',');
            value.forEach((primerStr) => {
              let queryObj = {
                Institution: primerStr.split('_')[0],
                Name: primerStr.split('_')[1],
              };
              const primer = queryPrimers(queryObj);
              if (primer !== undefined) {
                this.pendingChanges.configStore[param]
                  ? this.pendingChanges.configStore[param].push(primer)
                  : (this.pendingChanges.configStore[param] = [primer]);
              }
            });
          }
        } else if (param.toUpperCase() in GEO_LEVELS) {
          // Not in paramsToDisplay but can be set to/from URL
          // If a location is specified, update selectedLocationNodes
          const value = this.urlParams.getAll(param);

          value.forEach((item) => {
            const node = getLocationByNameAndLevel(
              rootStoreInstance.locationDataStore.selectTree,
              item,
              param,
              true
            )[0];

            if (typeof node !== 'undefined') {
              this.pendingChanges.configStore[param]
                ? this.pendingChanges.configStore[param].push(node)
                : (this.pendingChanges.configStore[param] = [node]);
            }
          });
        } else if (param === 'tab') {
          const value = this.urlParams.get(param);
          if (!Object.values(TABS).includes(value)) {
            // If not valid, set to home
            this.urlParams.set(param, TABS.TAB_EXAMPLE);
          }
        } else if (Object.keys(config['group_cols']).includes(param)) {
          // groupSelectFields
          const value = this.urlParams.get(param);
          this.pendingChanges.configStore[param]
            ? this.pendingChanges.configStore[param].push(value)
            : (this.pendingChanges.configStore[param] = [value]);
        } else {
          // A parameter was supplied in the url that cannot be set from the url
          // Remove that parameter from the url and do nothing
          this.urlParams.delete(param);
        }
      });
    });
    // Update any improperly formatted params
    updateURLFromParams(this.urlParams);
    // Update all stores
    Object.keys(this.pendingChanges).forEach((store) => {
      if (this.pendingChanges[store]) {
        this.updateStore(store);
      }
    });
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

  resetPendingChanges(store) {
    this.pendingChanges[store] = {};
  }

  updateStore(store) {
    // Propogate pendingChanges to store
    rootStoreInstance[store].applyPendingChanges(this.pendingChanges[store]);
    this.resetPendingChanges(store);
  }
}

import { action } from 'mobx';

import { rootStoreInstance } from './rootStore';

import { config } from '../config';
import {
  configStore as initialConfigStore,
  groupDataStore as initialGroupDataStore,
  plotSettingsStore as initialPlotSettingsStore,
} from '../constants/initialValues';

import { updateURLFromParams } from '../utils/updateQueryParam';
import { queryPrimers } from '../utils/primer';
import { getLocationByNameAndLevel } from '../utils/location';

import { textToCoords, textToResidueCoords } from '../utils/coordinates';

import { getGene, getProtein } from '../utils/gene_protein';

import { GEO_LEVELS, GROUP_MUTATION } from '../constants/defs.json';

export class URLMonitor {
  urlParams = new URLSearchParams(window.location.search);

  pendingChanges = {
    configStore: {},
    groupDataStore: {},
    plotSettingsStore: {},
  };

  init() {
    const defaultSelectedLocationNodes = [];
    if (initialConfigStore.defaultSelectedLocationNodes.length === 0) {
      initialConfigStore.defaultSelectedLocationNodes = ['USA', 'Canada'];
    }
    initialConfigStore.defaultSelectedLocationNodes.forEach((country) => {
      const countryNode = getLocationByNameAndLevel(
        rootStoreInstance.locationDataStore.selectTree,
        country,
        'country',
        true
      )[0];

      if (countryNode !== undefined) {
        defaultSelectedLocationNodes.push(countryNode);
      }
    });

    // FOR ANY ARRAY FIELDS -----
    // If the names of any array fields are in the URL params
    // then clear the array and set it to empty so we can add items later
    if (
      Object.keys(config.group_cols).some((groupKey) =>
        this.urlParams.has(groupKey)
      )
    ) {
      this.pendingChanges.configStore.selectedGroupFields = {};
    }
    if (this.urlParams.has('selectedPrimers')) {
      this.pendingChanges.configStore['selectedPrimers'] = [];
    }

    this.urlParams.forEach((value, key) => {
      value = decodeURIComponent(value);

      /* ------------
       * CONFIG STORE
       * ------------ */
      if (key in initialConfigStore) {
        if (key === 'selectedReference') {
          this.pendingChanges.configStore.selectedReference = value;
        }
        // Set gene/protein from URL params
        // Since the gene/protein object is dependent on the reference,
        // wait til we get that first and then we'll call getGene()/getProtein()
        // For now, set the value to a string
        else if (key === 'selectedGene') {
          const selectedReference = this.urlParams.has('selectedReference')
            ? this.urlParams.get('selectedReference')
            : rootStoreInstance.configStore.selectedReference;

          let selectedGene = getGene(key, selectedReference);
          if (selectedGene === undefined) {
            selectedGene = initialConfigStore.selectedGene;
          }

          this.pendingChanges.configStore[key] = selectedGene;
        } else if (key === 'selectedProtein') {
          const selectedReference = this.urlParams.has('selectedReference')
            ? this.urlParams.get('selectedReference')
            : rootStoreInstance.configStore.selectedReference;

          let selectedProtein = getProtein(key, selectedReference);
          if (selectedProtein === undefined) {
            selectedProtein = initialConfigStore.selectedProtein;
          }

          this.pendingChanges.configStore[key] = selectedProtein;
        } else if (key === 'ageRange' || key.includes('valid')) {
          // AgeRange is not being used currently so ignore
          // validity flags should not be set from the url
          return;
        } else if (key === 'customCoordinates') {
          this.pendingChanges.configStore[key] = textToCoords(value);
        } else if (key === 'residueCoordinates') {
          this.pendingChanges.configStore[key] = textToResidueCoords(value);
        } else if (key === 'customSequences') {
          // Store customSequences as an array of strings
          value = value.split(',');
          this.pendingChanges.configStore[key] = value;
        } else if (key === 'selectedPrimers') {
          value = value.split(',');
          value.forEach((primerStr) => {
            let queryObj = {
              Institution: primerStr.split('_')[0],
              Name: primerStr.split('_')[1],
            };
            const primer = queryPrimers(queryObj);
            if (primer !== undefined && primer !== null)
              this.pendingChanges.configStore[key].push(primer);
          });
        } else {
          this.pendingChanges.configStore[key] = value;
        }
      } else if (key in initialGroupDataStore) {
        /* ----------------
         * GROUP DATA STORE
         * ---------------- */
        if (key === 'selectedReportGroups') {
          value = decodeURIComponent(value).split(',');
          this.pendingChanges.groupDataStore[key] = value;
        } else {
          this.pendingChanges.groupDataStore[key] = value;
        }
      } else if (key in initialPlotSettingsStore) {
        /* -------------------
         * PLOT SETTINGS STORE
         * ------------------- */
        if (key === 'reportMutationListHidden') {
          value = decodeURIComponent(value).split(',');
          this.pendingChanges.plotSettingsStore[key] = value;
        } else {
          this.pendingChanges.plotSettingsStore[key] = value;
        }
      }

      // LOCATIONS
      else if (key.toUpperCase() in GEO_LEVELS) {
        if (
          !Object.prototype.hasOwnProperty.call(
            this.pendingChanges.configStore,
            'selectedLocationNodes'
          )
        ) {
          this.pendingChanges.configStore.selectedLocationNodes = [];
        }

        // If a location is specified, update selectedLocationNodes
        value = value.split(',');

        value.forEach((item) => {
          const node = getLocationByNameAndLevel(
            rootStoreInstance.locationDataStore.selectTree,
            item,
            key,
            true
          )[0];

          if (
            typeof node !== 'undefined' &&
            !this.pendingChanges.configStore.selectedLocationNodes.includes(
              node
            )
          ) {
            this.pendingChanges.configStore.selectedLocationNodes.push(node);
          }
        });
      }
      // selectedGroupFields
      else if (Object.keys(config['group_cols']).includes(key)) {
        if (
          !Object.prototype.hasOwnProperty.call(
            this.pendingChanges.configStore,
            'selectedGroupFields'
          )
        ) {
          this.pendingChanges.configStore.selectedGroupFields = {};
        }

        if (
          !Object.prototype.hasOwnProperty.call(
            this.pendingChanges.configStore.selectedGroupFields,
            key
          )
        ) {
          this.pendingChanges.configStore.selectedGroupFields[key] = [];
        }
        if (
          !this.pendingChanges.configStore.selectedGroupFields[key].includes(
            value
          )
        ) {
          this.pendingChanges.configStore.selectedGroupFields[key].push(value);
        }
      } else if (key === 'subtype') {
        if (
          !Object.prototype.hasOwnProperty.call(
            this.pendingChanges.configStore,
            'selectedGroupFields'
          )
        ) {
          this.pendingChanges.configStore.selectedGroupFields = {};
        }
        this.pendingChanges.configStore.selectedGroupFields[key] = [value];
      }
      // Ignore for now -- only change tab after all stores have
      // pending changes from URL params flushed to their state
      // For now - just don't delete this param from the URL
      else if (key === 'tab') {
        return;
      } else {
        // Invalid field, remove it from the url
        console.warn('DELETING URL PARAM ' + key);
        this.urlParams.delete(key);
      }
    });

    // If no locations in url, set default selected locations
    if (
      !Object.prototype.hasOwnProperty.call(
        this.pendingChanges.configStore,
        'selectedLocationNodes'
      ) ||
      this.pendingChanges.configStore.selectedLocationNodes.length === 0
    ) {
      this.pendingChanges.configStore.selectedLocationNodes =
        defaultSelectedLocationNodes;
    }

    // RSV EDGE CASE
    // If we don't start out in mutation mode, and no
    // selected group fields are specified, empty out the
    // selected group fields selection
    const groupKey = this.urlParams.has('groupKey')
      ? this.urlParams.get('groupKey')
      : rootStoreInstance.configStore.groupKey;
    if (
      groupKey !== GROUP_MUTATION &&
      Object.keys(config.group_cols).every((key) => !this.urlParams.has(key))
    ) {
      this.pending.configStore.selectedGroupFields = {};
    }

    // Update URL
    updateURLFromParams(this.urlParams);

    // Update all stores
    Object.keys(this.pendingChanges).forEach((store) => {
      this.updateStore(store);
    });

    // HANDLE TAB SELECTION
    if (this.urlParams.has('tab')) {
      rootStoreInstance.UIStore.setActiveTab(this.urlParams.get('tab'));
    }
  }

  @action
  updateStore(store) {
    // console.log('updateStore', store, this.pendingChanges[store]);
    // Propogate pendingChanges to store
    rootStoreInstance[store].applyPendingChanges(
      this.pendingChanges[store],
      false
    );
    rootStoreInstance[store].updateURL(this.pendingChanges[store]);
    // Reset pending changes
    this.pendingChanges[store] = {};
  }
}

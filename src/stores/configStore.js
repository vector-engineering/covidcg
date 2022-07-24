import {
  observable,
  action,
  toJS,
  //intercept, autorun
} from 'mobx';

import { getGene, getProtein } from '../utils/gene_protein';
import { queryReferenceSequence } from '../utils/reference';
import { getLocationByNameAndLevel } from '../utils/location';
import { updateURLFromParams } from '../utils/updateQueryParam';
import { queryPrimers } from '../utils/primer';
import { arrayEqual } from '../utils/func';

import {
  GROUP_MUTATION,
  DNA_OR_AA,
  COORDINATE_MODES,
  GEO_LEVELS,
  TABS,
  GROUPS,
} from '../constants/defs.json';
import { config } from '../config';

import { PARAMS_TO_TRACK } from './paramsToTrack';
import { rootStoreInstance } from './rootStore';
import { configStore as initialConfigStore } from '../constants/initialValues';

const urlParams = new URLSearchParams(window.location.search);

const defaultsFromParams = {};

PARAMS_TO_TRACK.forEach((param) => {
  // console.log('getting: ', param, urlParams.get(param));
  defaultsFromParams[param] = urlParams.get(param);
});

export class ConfigStore {
  // Maintain a reference to the initial values
  // Initalize values to correct data types
  initialValues = {};

  @observable groupKey = '';
  @observable dnaOrAa = '';

  @observable selectedReference = '';

  @observable selectedGene = {};
  @observable selectedProtein = {};
  @observable selectedPrimers = [];

  @observable customCoordinates = [[]];
  @observable customSequences = [];
  @observable residueCoordinates = [[]];
  @observable coordinateMode = '';

  @observable startDate = new Date();
  @observable endDate = new Date();

  @observable submStartDate = '';
  @observable submEndDate = '';

  @observable selectedGroupFields = {};
  @observable selectedLocationNodes = [];

  @observable hoverGroup = null;
  @observable selectedGroups = [];

  @observable selectedMetadataFields = {};
  @observable ageRange = [];

  @observable hoverLocation = null;
  @observable focusedLocations = [];

  constructor() {}

  init() {
    this.initialValues = initialConfigStore;

    Object.keys(this.initialValues).forEach((key) => {
      this[key] = this.initialValues[key];
    });

    PARAMS_TO_TRACK.forEach((param) => {
      if (defaultsFromParams[param]) {
        // console.log('setting: ', param, urlParams.get(param));
        this[param] = defaultsFromParams[param];
      }
    });

    this.urlParams = new URLSearchParams(window.location.search);

    // RSV EDGE CASE
    // If any selected group fields are set in the URLs,
    // then clear the default selected group fields
    if (
      Object.keys(config.group_cols).some((groupKey) =>
        this.urlParams.has(groupKey)
      )
    ) {
      this.selectedGroupFields = {};
    }

    // Check to see what's in the URL
    this.urlParams.forEach((value, key) => {
      value = decodeURIComponent(value);
      if (key in this.initialValues) {
        if (key === 'selectedReference') {
          this.selectedReference = value;
        }
        // Set gene/protein from URL params
        // Since the gene/protein object is dependent on the reference,
        // wait til we get that first and then we'll call getGene()/getProtein()
        // For now, set the value to a string
        else if (key === 'selectedGene' || key === 'selectedProtein') {
          this[key] = value;
        } else if (key === 'ageRange' || key.includes('valid')) {
          // AgeRange is not being used currently so ignore
          // validity flags should not be set from the url
          return;
        } else if (
          key === 'customCoordinates' ||
          key === 'residueCoordinates'
        ) {
          // If coordinates are specified, save them as an array of numbers
          // Coordinates can stay as a string in the URL
          let arr = [];
          value = value.split(',');
          value.forEach((item, i) => {
            value[i] = parseInt(item);

            if (i % 2 === 1) {
              arr.push([value[i - 1], value[i]]);
            }
          });

          this[key] = arr;
        } else if (key === 'customSequences') {
          // Store customSequences as an array of strings
          value = value.split(',');
          this[key] = value;
        } else if (key === 'selectedPrimers') {
          value = value.split(',');
          value.forEach((primerStr) => {
            let queryObj = {
              Institution: primerStr.split('_')[0],
              Name: primerStr.split('_')[1],
            };
            const primer = queryPrimers(queryObj);
            if (primer !== undefined) this[key].push(primer);
          });
        } else {
          this[key] = value;
        }
      } else if (key.toUpperCase() in GEO_LEVELS) {
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
            !this.selectedLocationNodes.includes(node)
          ) {
            this.selectedLocationNodes.push(node);
          }
        });
      } else if (key === 'tab') {
        // Check if the specified tab value is valid (included in TABS)
        // tab is read and activeTab is set from routes.js
        if (Object.values(TABS).includes(value)) {
          this[key] = value;
        } else {
          // If not valid, set to home
          this.urlParams.set(key, TABS.TAB_EXAMPLE);
          this[key] = TABS.TAB_EXAMPLE;
        }
      }
      // selectedGroupFields
      else if (Object.keys(config['group_cols']).includes(key)) {
        if (
          !Object.prototype.hasOwnProperty.call(this.selectedGroupFields, key)
        ) {
          this.selectedGroupFields[key] = [];
        }
        if (!this.selectedGroupFields[key].includes(value)) {
          this.selectedGroupFields[key].push(value);
        }
      } else if (key === 'subtype') {
        this.selectedGroupFields[key] = [value];
      } else {
        // Invalid field, remove it from the url
        // console.log('DELETE ' + key);
        this.urlParams.delete(key);
      }
    });

    // Map gene/protein names to objects
    if (typeof this.selectedGene === 'string') {
      this.selectedGene = getGene(this.selectedGene, this.selectedReference);
      if (this.selectedGene === undefined) {
        this.selectedGene = this.initialValues.selectedGene;
      }
    }
    if (typeof this.selectedProtein === 'string') {
      this.selectedProtein = getProtein(
        this.selectedProtein,
        this.selectedReference
      );
      if (this.selectedProtein === undefined) {
        this.selectedProtein = this.initialValues.selectedProtein;
      }
    }

    // Update URL
    updateURLFromParams(this.urlParams);

    const defaultSelectedLocationNodes = [
      getLocationByNameAndLevel(
        rootStoreInstance.locationDataStore.selectTree,
        'USA',
        'country',
        true
      )[0],
      getLocationByNameAndLevel(
        rootStoreInstance.locationDataStore.selectTree,
        'Canada',
        'country',
        true
      )[0],
    ].filter((node) => node !== undefined);
    this.initialValues['selectedLocationNodes'] = defaultSelectedLocationNodes;

    if (this.selectedLocationNodes.length === 0) {
      // If no locations in url, set default selected locations
      this.selectedLocationNodes = defaultSelectedLocationNodes;
    }

    // RSV EDGE CASE
    // If we don't start out in mutation mode, and no
    // selected group fields are specified, empty out the
    // selected group fields selection
    if (
      this.groupKey !== GROUP_MUTATION &&
      Object.keys(config.group_cols).every(
        (groupKey) => !this.urlParams.has(groupKey)
      )
    ) {
      this.selectedGroupFields = {};
    }
  }

  // modifyQueryParams = autorun(() => {
  //   PARAMS_TO_TRACK.forEach((param) => {
  //     updateQueryStringParam(param, JSON.stringify(this[param]));
  //   });
  // });

  @action
  resetValues = (values) => {
    Object.keys(this.initialValues).forEach((key) => {
      if (key in values) {
        this[key] = values[key];
      } else {
        this[key] = this.initialValues[key];
      }

      // Special actions for some keys
      if (key === 'selectedLocationNodes') {
        rootStoreInstance.locationDataStore.setSelectedNodes(values[key]);
      }
    });

    // Manually set residue coordinates if they weren't specified
    if (
      'selectedGene' in values &&
      !('residueCoordinates' in values) &&
      this.selectedGene.name !== 'All Genes'
    ) {
      this.residueCoordinates = [[1, this.selectedGene.len_aa]];
    } else if (
      'selectedProtein' in values &&
      !('residueCoordinates' in values) &&
      this.selectedProtein.name !== 'All Proteins'
    ) {
      this.residueCoordinates = [[1, this.selectedProtein.len_aa]];
    }

    // Trigger data re-run
    rootStoreInstance.dataStore.fetchData();
  };

  @action
  applyPendingChanges = (pending) => {
    // Clear selected groups/locations
    this.hoverGroup = this.initialValues.hoverGroup;
    this.selectedGroups = this.initialValues.selectedGroups;
    this.hoverLocation = this.initialValues.hoverLocation;
    this.focusedLocations = this.initialValues.focusedLocations;

    // Overwrite any of our fields here with the pending ones
    Object.keys(pending).forEach((field) => {
      this[field] = pending[field];

      // Update urlParams
      this.urlParams.delete(field);

      if (field === 'selectedGene' || field === 'selectedProtein') {
        // Handle fields that return objects
        this.urlParams.set(field, pending[field].name);
      } else if (
        field === 'selectedMetadataFields' ||
        field === 'selectedLocationNodes' ||
        field === 'ageRange'
      ) {
        // Ignore Metadata fields
        // selectedLocationNodes is displayed as node.level=node.value in the URL
        // ageRange is not currently being used
        return;
      } else if (field.includes('valid')) {
        // Ignore boolean flags
        return;
      } else if (field === 'selectedPrimers') {
        pending[field].forEach((primer) => {
          if (this.urlParams.has(field)) {
            this.urlParams.append(
              field,
              primer.Institution + '_' + primer.Name
            );
          } else {
            this.urlParams.set(field, primer.Institution + '_' + primer.Name);
          }
        });
      } else if (field === 'selectedGroupFields') {
        // If selectedGroupFields just got emptied - make sure to
        // remove any lingering group keys in the URL
        Object.keys(config.group_cols).forEach((groupKey) => {
          if (!Object.keys(pending[field]).includes(groupKey)) {
            this.urlParams.delete(groupKey);
          }
        });

        Object.keys(pending[field]).forEach((groupKey) => {
          this.urlParams.delete(groupKey);
          pending[field][groupKey].forEach((group) => {
            this.urlParams.append(groupKey, group);
          });
        });
      } else {
        this.urlParams.set(field, String(pending[field]));
      }

      if (pending[field] === this.initialValues[field]) {
        // Only display non-default fields in the url
        this.urlParams.delete(field);
      }
    });

    // Show only relevant coordinate info
    const mode = this.urlParams.get('coordinateMode');
    if (mode === 'protein') {
      this.urlParams.delete('selectedGene');
      this.urlParams.delete('selectedPrimers');
      this.urlParams.delete('customCoordinates');
      this.urlParams.delete('customSequences');
    } else if (mode === 'gene' || !mode) {
      // Gene is the default mode and may have been deleted so check for null
      this.urlParams.delete('selectedProtein');
      this.urlParams.delete('selectedPrimers');
      this.urlParams.delete('customCoordinates');
      this.urlParams.delete('customSequences');
    } else if (mode === 'primer') {
      this.urlParams.delete('selectedProtein');
      this.urlParams.delete('selectedGene');
      this.urlParams.delete('customCoordinates');
      this.urlParams.delete('customSequences');
    } else if (mode === 'custom') {
      this.urlParams.delete('selectedProtein');
      this.urlParams.delete('selectedPrimers');
      this.urlParams.delete('selectedGene');
      this.urlParams.delete('customSequences');
    } else if (mode === 'sequence') {
      this.urlParams.delete('selectedProtein');
      this.urlParams.delete('selectedPrimers');
      this.urlParams.delete('selectedGene');
      this.urlParams.delete('customCoordinates');
    }

    if (this.groupKey !== GROUP_MUTATION) {
      this.urlParams.delete('residueCoordinates');
      this.urlParams.delete('selectedGene');
      this.urlParams.delete('selectedProtein');
      this.urlParams.delete('selectedPrimers');
      this.urlParams.delete('customCoordinates');
      this.urlParams.delete('customSequences');
    }

    // Update the location node tree with our new selection
    rootStoreInstance.locationDataStore.setSelectedNodes(
      this.selectedLocationNodes
    );

    // Update the location URL params
    Object.values(GEO_LEVELS).forEach((level) => {
      this.urlParams.delete(level);
    });

    this.selectedLocationNodes.forEach((node) => {
      if (this.urlParams.has(String(node.level))) {
        this.urlParams.append(String(node.level), String(node.value_txt));
      } else {
        this.urlParams.set(String(node.level), String(node.value_txt));
      }
    });

    updateURLFromParams(this.urlParams);

    // Get the new data from the server
    rootStoreInstance.dataStore.fetchData();
  };

  getMutationType() {
    if (this.dnaOrAa === DNA_OR_AA.DNA) {
      return 'dna';
    } else if (this.dnaOrAa === DNA_OR_AA.AA) {
      if (this.coordinateMode === COORDINATE_MODES.COORD_GENE) {
        return 'gene_aa';
      } else if (this.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
        return 'protein_aa';
      }
    }
    return undefined;
  }

  // Get a pretty name for the group
  getGroupLabel(groupKey = null, dnaOrAa = null) {
    // Default to using store attribute, if no explicit groupKey
    // is provided
    if (groupKey === null) {
      groupKey = this.groupKey;
    }
    if (dnaOrAa === null) {
      dnaOrAa = this.dnaOrAa;
    }

    if (Object.keys(config.group_cols).includes(groupKey)) {
      return config.group_cols[groupKey].title;
    } else if (groupKey === GROUP_MUTATION) {
      if (dnaOrAa === DNA_OR_AA.DNA) {
        return 'NT Mutation';
      } else {
        return 'AA Mutation';
      }
    }
  }

  getSelectedLocations() {
    const res = {
      region: [],
      country: [],
      division: [],
      location: [],
    };
    //console.log(this.selectedLocationNodes);
    this.selectedLocationNodes.forEach((node) => {
      res[node.level].push(node.value);
    });
    return res;
  }

  getCoordinateRanges() {
    // Set the coordinate range based off the coordinate mode
    if (this.coordinateMode === COORDINATE_MODES.COORD_GENE) {
      // Return ranges if All Genes
      if (this.selectedGene.name === 'All Genes') {
        return this.selectedGene.ranges;
      }
      // Disable residue indices for non-protein-coding genes
      if (!this.selectedGene.protein_coding) {
        return this.selectedGene.segments;
      }
      const coordinateRanges = [];
      this.residueCoordinates.forEach((range) => {
        // Make a deep copy of the current range
        const curRange = range.slice();

        if (this.dnaOrAa === DNA_OR_AA.DNA) {
          for (let i = 0; i < this.selectedGene.aa_ranges.length; i++) {
            const curAARange = this.selectedGene.aa_ranges[i];
            const curNTRange = this.selectedGene.segments[i];
            if (
              (curRange[0] >= curAARange[0] && curRange[0] <= curAARange[1]) ||
              (curRange[0] <= curAARange[0] && curRange[1] >= curAARange[0])
            ) {
              coordinateRanges.push([
                curNTRange[0] + (curRange[0] - curAARange[0]) * 3,
                curNTRange[0] -
                  1 +
                  Math.min(curRange[1] - curAARange[0] + 1, curAARange[1]) * 3,
              ]);
              // Push the beginning of the current range to the end of
              // the current AA range of the gene
              if (curAARange[1] < curRange[1]) {
                curRange[0] = curAARange[1] + 1;
              }
            }
          }
        } else {
          coordinateRanges.push([curRange[0], curRange[1]]);
        }
      });
      return coordinateRanges;
    } else if (this.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
      const coordinateRanges = [];
      this.residueCoordinates.forEach((range) => {
        // Make a deep copy of the current range
        const curRange = range.slice();

        if (this.dnaOrAa === DNA_OR_AA.DNA) {
          for (let i = 0; i < this.selectedProtein.aa_ranges.length; i++) {
            const curAARange = this.selectedProtein.aa_ranges[i];
            const curNTRange = this.selectedProtein.segments[i];
            if (
              (curRange[0] >= curAARange[0] && curRange[0] <= curAARange[1]) ||
              (curRange[0] <= curAARange[0] && curRange[1] >= curAARange[0])
            ) {
              coordinateRanges.push([
                curNTRange[0] + (curRange[0] - curAARange[0]) * 3,
                curNTRange[0] -
                  1 +
                  Math.min(curRange[1] - curAARange[0] + 1, curAARange[1]) * 3,
              ]);
              // Push the beginning of the current range to the end of
              // the current AA range of the gene
              if (curAARange[1] < curRange[1]) {
                curRange[0] = curAARange[1] + 1;
              }
            }
          }
        } else {
          coordinateRanges.push([curRange[0], curRange[1]]);
        }
      });
      return coordinateRanges;
    } else if (this.coordinateMode === COORDINATE_MODES.COORD_PRIMER) {
      return this.selectedPrimers.map((primer) => {
        return [primer.Start, primer.End];
      });
    } else if (this.coordinateMode === COORDINATE_MODES.COORD_CUSTOM) {
      return toJS(this.customCoordinates);
    } else if (this.coordinateMode === COORDINATE_MODES.COORD_SEQUENCE) {
      return this.customSequences.map((seq) => {
        return queryReferenceSequence(seq, this.selectedReference);
      });
    }
  }

  getSelectedMetadataFields() {
    const selectedMetadataFields = toJS(this.selectedMetadataFields);
    Object.keys(selectedMetadataFields).forEach((metadataField) => {
      selectedMetadataFields[metadataField] = selectedMetadataFields[
        metadataField
      ].map((item) => {
        return parseInt(item.value);
      });
    });
    return selectedMetadataFields;
  }

  @action
  updateHoverGroup = (group) => {
    // console.log('UPDATE HOVER GROUP', group);
    if (group === this.hoverGroup) {
      return;
    } else if (group === GROUPS.NONE_GROUP || GROUPS.ALL_OTHER_GROUP) {
      // Ignore for some special groups
      return;
    } else if (group === null) {
      this.hoverGroup = null;
    } else {
      this.hoverGroup = group;
    }
  };

  @action
  updateSelectedGroups = (groups) => {
    // First check to see that it's different. If not,
    // skip the update
    // This matters because JS arrays are passed by reference,
    // and any listeners which see a change to this reference
    // will update themselves when the reference changes,
    // even if the underlying data does not
    // groups are in the structure [{ group: group }]
    if (
      arrayEqual(
        groups.map((group) => group.group),
        this.selectedGroups.map((group) => group.group)
      )
    ) {
      return;
    }

    this.selectedGroups = groups;
    if (this.groupKey === GROUP_MUTATION) {
      rootStoreInstance.dataStore.processSelectedMutations();
    }
  };

  getSelectedGroupIds() {
    const { dnaMutationMap, geneAaMutationMap, proteinAaMutationMap } =
      rootStoreInstance.mutationDataStore;

    let selectedGroupIds;
    if (this.dnaOrAa === DNA_OR_AA.DNA) {
      selectedGroupIds = this.selectedGroups
        .map((item) => dnaMutationMap[item.group])
        .map((mutationId) =>
          mutationId === undefined ? -1 : parseInt(mutationId)
        );
    } else if (this.dnaOrAa === DNA_OR_AA.AA) {
      if (this.coordinateMode === COORDINATE_MODES.COORD_GENE) {
        selectedGroupIds = this.selectedGroups
          .map((item) => geneAaMutationMap[item.group])
          .map((mutationId) =>
            mutationId === undefined ? -1 : parseInt(mutationId)
          );
      } else if (this.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
        selectedGroupIds = this.selectedGroups
          .map((item) => proteinAaMutationMap[item.group])
          .map((mutationId) =>
            mutationId === undefined ? -1 : parseInt(mutationId)
          );
      }
    }
    // Array to Set
    selectedGroupIds = new Set(selectedGroupIds);

    return selectedGroupIds;
  }

  getIntToMutationMap() {
    const {
      intToDnaMutationMap,
      intToGeneAaMutationMap,
      intToProteinAaMutationMap,
    } = rootStoreInstance.mutationDataStore;

    if (this.dnaOrAa === DNA_OR_AA.DNA) {
      return intToDnaMutationMap;
    } else if (this.dnaOrAa === DNA_OR_AA.AA) {
      if (this.coordinateMode === COORDINATE_MODES.COORD_GENE) {
        return intToGeneAaMutationMap;
      } else if (this.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
        return intToProteinAaMutationMap;
      }
    }
  }

  getMutationToIntMap() {
    const { dnaMutationMap, geneAaMutationMap, proteinAaMutationMap } =
      rootStoreInstance.mutationDataStore;

    if (this.dnaOrAa === DNA_OR_AA.DNA) {
      return dnaMutationMap;
    } else if (this.dnaOrAa === DNA_OR_AA.AA) {
      if (this.coordinateMode === COORDINATE_MODES.COORD_GENE) {
        return geneAaMutationMap;
      } else if (this.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
        return proteinAaMutationMap;
      }
    }
  }

  @action
  updateHoverLocation = (location) => {
    this.hoverLocation = location;
  };

  @action
  updateFocusedLocations = (locations) => {
    if (
      arrayEqual(
        locations.map((location) => location.location),
        this.focusedLocations.map((location) => location.location)
      )
    ) {
      return;
    }
    this.focusedLocations = locations;
  };
}

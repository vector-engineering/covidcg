import {
  observable,
  action,
  toJS,
  //intercept, autorun
} from 'mobx';

import {
  geneMap,
  proteinMap,
  getGene,
  getProtein,
} from '../utils/gene_protein';
import { queryReferenceSequence } from '../utils/reference';
import { getLocationByNameAndLevel } from '../utils/location';
import { intToISO, ISOToInt } from '../utils/date';
import { updateURLFromParams } from '../utils/updateQueryParam';
import { queryPrimers } from '../utils/primer';
import { arrayEqual } from '../utils/func';

import {
  GROUP_SNV,
  DNA_OR_AA,
  COORDINATE_MODES,
  GEO_LEVELS,
  TABS,
  GROUPS,
} from '../constants/defs.json';
import { config } from '../config';

// import { updateQueryStringParam } from '../utils/updateQueryParam';
import { PARAMS_TO_TRACK } from './paramsToTrack';
import { rootStoreInstance } from './rootStore';

// Define initial values

const today = intToISO(new Date().getTime());
const lastNDays = 30; // By default, show only the last 1 month

export const initialValues = {
  groupKey: 'snv',
  dnaOrAa: DNA_OR_AA.AA,

  // Select the Spike gene and nsp13 protein by default
  selectedGene: getGene('S'),
  selectedProtein: getProtein('nsp12 - RdRp'),
  selectedPrimers: [],
  customCoordinates: [[8000, 12000]],
  customSequences: ['GACCCCAAAATCAGCGAAAT'],
  residueCoordinates: [[1, getGene('S').len_aa]],

  // Selecting the gene as the coordinate range by default
  coordinateMode: COORDINATE_MODES.COORD_GENE,

  // days * (24 hours/day) * (60 min/hour) * (60 s/min) * (1000 ms/s)
  startDate: intToISO(ISOToInt(today) - lastNDays * 24 * 60 * 60 * 1000),
  endDate: today,

  submStartDate: '',
  submEndDate: '',

  selectedLocationNodes: [],

  hoverGroup: null,
  selectedGroups: [],

  // Metadata filtering
  selectedMetadataFields: {},
  ageRange: [null, null],

  // Location tab
  hoverLocation: null,
  focusedLocations: [],
};

const urlParams = new URLSearchParams(window.location.search);

const defaultsFromParams = {};

PARAMS_TO_TRACK.forEach((param) => {
  // console.log('getting: ', param, urlParams.get(param));
  defaultsFromParams[param] = urlParams.get(param);
});

export class ConfigStore {
  // Maintain a reference to the initial values
  initialValues = initialValues;

  @observable groupKey = initialValues.groupKey;
  @observable dnaOrAa = initialValues.dnaOrAa;

  @observable selectedGene = initialValues.selectedGene;
  @observable selectedProtein = initialValues.selectedProtein;
  @observable selectedPrimers = initialValues.selectedPrimers;

  @observable customCoordinates = initialValues.customCoordinates;
  @observable customSequences = initialValues.customSequences;
  @observable residueCoordinates = initialValues.residueCoordinates;
  @observable coordinateMode = initialValues.coordinateMode;

  @observable startDate = initialValues.startDate;
  @observable endDate = initialValues.endDate;

  @observable submStartDate = initialValues.submStartDate;
  @observable submEndDate = initialValues.submEndDate;

  @observable selectedLocationNodes = initialValues.selectedLocationNodes;

  @observable hoverGroup = initialValues.hoverGroup;
  @observable selectedGroups = initialValues.selectedGroups;

  @observable selectedMetadataFields = initialValues.selectedMetadataFields;
  @observable ageRange = initialValues.ageRange;

  @observable hoverLocation = initialValues.hoverLocation;
  @observable focusedLocations = initialValues.focusedLocations;

  constructor() {}

  init() {
    PARAMS_TO_TRACK.forEach((param) => {
      if (defaultsFromParams[param]) {
        // console.log('setting: ', param, urlParams.get(param));
        this[param] = defaultsFromParams[param];
      }
    });

    this.urlParams = new URLSearchParams(window.location.search);

    // Check to see what's in the URL
    this.urlParams.forEach((value, key) => {
      value = decodeURIComponent(value);
      if (key in initialValues) {
        if (key === 'selectedGene') {
          // If the specified gene is in the geneMap get the gene
          if (value in geneMap) {
            this[key] = getGene(value);
          } else {
            // Else display default gene
            this[key] = initialValues.selectedGene;
            this.urlParams.set(key, initialValues.selectedGene.name);
          }
        } else if (key === 'selectedProtein') {
          // If the specified protein is in the proteinMap get the protein
          if (value in proteinMap) {
            this[key] = getProtein(value);
          } else {
            // Else display default protein
            this[key] = initialValues.selectedProtein;
            this.urlParams.set(key, initialValues.selectedProtein.name);
          }
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
      } else {
        // Invalid field, remove it from the url
        this.urlParams.delete(key);
      }
    });

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
  }

  // modifyQueryParams = autorun(() => {
  //   PARAMS_TO_TRACK.forEach((param) => {
  //     updateQueryStringParam(param, JSON.stringify(this[param]));
  //   });
  // });

  @action
  resetValues = (values) => {
    Object.keys(initialValues).forEach((key) => {
      if (key in values) {
        this[key] = values[key];
      } else {
        this[key] = initialValues[key];
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
    this.hoverGroup = initialValues.hoverGroup;
    this.selectedGroups = initialValues.selectedGroups;
    this.hoverLocation = initialValues.hoverLocation;
    this.focusedLocations = initialValues.focusedLocations;

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
      } else {
        this.urlParams.set(field, String(pending[field]));
      }

      if (pending[field] === initialValues[field]) {
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

  getSnvType() {
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
    } else if (groupKey === GROUP_SNV) {
      if (dnaOrAa === DNA_OR_AA.DNA) {
        return 'NT SNV';
      } else {
        return 'AA SNV';
      }
    }
  }

  getSelectedLocations() {
    const res = {
      region: [],
      country: [],
      division: [],
      location: []
    };
    this.selectedLocationNodes.forEach((node) => {
      res[node.level].push(node.value);
    });
    return res
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
      });
      return coordinateRanges;
    } else if (this.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
      const coordinateRanges = [];
      this.residueCoordinates.forEach((range) => {
        // Make a deep copy of the current range
        const curRange = range.slice();
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
        return queryReferenceSequence(seq);
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
    if (this.groupKey === GROUP_SNV) {
      rootStoreInstance.dataStore.processSelectedSnvs();
    }
  };

  getSelectedGroupIds() {
    const {
      dnaSnvMap,
      geneAaSnvMap,
      proteinAaSnvMap,
    } = rootStoreInstance.snpDataStore;

    let selectedGroupIds;
    if (this.dnaOrAa === DNA_OR_AA.DNA) {
      selectedGroupIds = this.selectedGroups
        .map((item) => dnaSnvMap[item.group])
        .map((snpId) => (snpId === undefined ? -1 : parseInt(snpId)));
    } else if (this.dnaOrAa === DNA_OR_AA.AA) {
      if (this.coordinateMode === COORDINATE_MODES.COORD_GENE) {
        selectedGroupIds = this.selectedGroups
          .map((item) => geneAaSnvMap[item.group])
          .map((snpId) => (snpId === undefined ? -1 : parseInt(snpId)));
      } else if (this.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
        selectedGroupIds = this.selectedGroups
          .map((item) => proteinAaSnvMap[item.group])
          .map((snpId) => (snpId === undefined ? -1 : parseInt(snpId)));
      }
    }
    // Array to Set
    selectedGroupIds = new Set(selectedGroupIds);

    return selectedGroupIds;
  }

  getIntToSnvMap() {
    const {
      intToDnaSnvMap,
      intToGeneAaSnvMap,
      intToProteinAaSnvMap,
    } = rootStoreInstance.snpDataStore;

    if (this.dnaOrAa === DNA_OR_AA.DNA) {
      return intToDnaSnvMap;
    } else if (this.dnaOrAa === DNA_OR_AA.AA) {
      if (this.coordinateMode === COORDINATE_MODES.COORD_GENE) {
        return intToGeneAaSnvMap;
      } else if (this.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
        return intToProteinAaSnvMap;
      }
    }
  }

  getSnvToIntMap() {
    const {
      dnaSnvMap,
      geneAaSnvMap,
      proteinAaSnvMap,
    } = rootStoreInstance.snpDataStore;

    if (this.dnaOrAa === DNA_OR_AA.DNA) {
      return dnaSnvMap;
    } else if (this.dnaOrAa === DNA_OR_AA.AA) {
      if (this.coordinateMode === COORDINATE_MODES.COORD_GENE) {
        return geneAaSnvMap;
      } else if (this.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
        return proteinAaSnvMap;
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

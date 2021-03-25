import {
  observable,
  action,
  toJS,
  //intercept, autorun
} from 'mobx';

import { geneMap, proteinMap,
         getGene, getProtein } from '../utils/gene_protein';
import { queryReferenceSequence } from '../utils/reference';
import { getLocationByNameAndLevel } from '../utils/location';
import { intToISO } from '../utils/date';
import { updateURLFromParams } from '../utils/parseQueryParams';

import {
  GROUP_SNV,
  DNA_OR_AA,
  COORDINATE_MODES,
  LOW_FREQ_FILTER_TYPES,
  MIN_DATE,
  COLOR_MODES,
  COMPARE_MODES,
  COMPARE_COLORS,
  TABS,
  GEO_LEVELS,
} from '../constants/defs.json';
import { config } from '../config';

// import { updateQueryStringParam } from '../utils/updateQueryParam';
import { PARAMS_TO_TRACK } from './paramsToTrack';
import { rootStoreInstance } from './rootStore';

// Define initial values

export const initialConfigValues = {
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

  startDate: MIN_DATE,
  endDate: intToISO(new Date().getTime()),

  selectedLocationNodes: [],

  hoverGroup: null,
  selectedGroups: [],

  // Metadata filtering
  selectedMetadataFields: {},
  ageRange: [null, null],

  // Location tab
  hoverLocation: null,
  focusedLocations: [],

  lowFreqFilterType: LOW_FREQ_FILTER_TYPES.GROUP_COUNTS,
  maxGroupCounts: 100,
  minLocalCounts: 50,
  minGlobalCounts: 100,

  // Query String
  urlParams: {
    tab: 'home',
  },
};

const urlParams = new URLSearchParams(window.location.search);

const defaultsFromParams = {};

PARAMS_TO_TRACK.forEach((param) => {
  // console.log('getting: ', param, urlParams.get(param));
  defaultsFromParams[param] = urlParams.get(param);
});

export class ConfigStore {
  // Maintain a reference to the initial values
  initialConfigValues = initialConfigValues;

  // References to store instances
  plotSettingsStoreInstance;
  locationDataStoreInstance;
  dataStoreInstance;

  @observable groupKey = initialConfigValues.groupKey;
  @observable dnaOrAa = initialConfigValues.dnaOrAa;

  @observable selectedGene = initialConfigValues.selectedGene;
  @observable selectedProtein = initialConfigValues.selectedProtein;
  @observable selectedPrimers = initialConfigValues.selectedPrimers;

  @observable customCoordinates = initialConfigValues.customCoordinates;
  @observable customSequences = initialConfigValues.customSequences;
  @observable residueCoordinates = initialConfigValues.residueCoordinates;
  @observable coordinateMode = initialConfigValues.coordinateMode;

  @observable startDate = initialConfigValues.startDate;
  @observable endDate = initialConfigValues.endDate;

  @observable selectedLocationNodes = initialConfigValues.selectedLocationNodes;

  @observable hoverGroup = initialConfigValues.hoverGroup;
  @observable selectedGroups = initialConfigValues.selectedGroups;

  @observable selectedMetadataFields =
    initialConfigValues.selectedMetadataFields;
  @observable ageRange = initialConfigValues.ageRange;

  @observable hoverLocation = initialConfigValues.hoverLocation;
  @observable focusedLocations = initialConfigValues.focusedLocations;

  @observable lowFreqFilterType = initialConfigValues.lowFreqFilterType;
  @observable maxGroupCounts = initialConfigValues.maxGroupCounts;
  @observable minLocalCounts = initialConfigValues.minLocalCounts;
  @observable minGlobalCounts = initialConfigValues.minGlobalCounts;

  urlParams = initialConfigValues.urlParams;

  constructor() {}

  init() {
    this.plotSettingsStoreInstance = rootStoreInstance.plotSettingsStore;
    this.locationDataStoreInstance = rootStoreInstance.locationDataStore;
    this.dataStoreInstance = rootStoreInstance.dataStore;
    this.snpDataStoreInstance = rootStoreInstance.snpDataStore;

    PARAMS_TO_TRACK.forEach((param) => {
      if (defaultsFromParams[param]) {
        // console.log('setting: ', param, urlParams.get(param));
        this[param] = defaultsFromParams[param];
      }
    });

    this.urlParams = new URLSearchParams(window.location.search);

    // Check to see what's in the URL
    this.urlParams.forEach((value, key) => {
      if (key.toUpperCase() in GEO_LEVELS) {
        // If a location is specified, update selectedLocationNodes
        value = encodeURIComponent(value);
        value = value.split('%2C');

        value.forEach((item) => {
          const node = getLocationByNameAndLevel(
            this.locationDataStoreInstance.selectTree,
            item,
            key,
            true
          )[0];

          if (typeof node !== 'undefined') {
            this.selectedLocationNodes.push(node);
          }
        });
      } else if (key === 'selectedGene') {
        // If the specified gene is in the geneMap get the gene
        if (value in geneMap) {
          this[key] = getGene(value);
        } else {
          // Else display default gene
          this[key] = getGene('S');
          this.urlParams.set(key, 'S');
        }
      } else if (key === 'selectedProtein') {
        // If the specified protein is in the proteinMap get the protein
        if (value in proteinMap) {
          this[key] = getProtein(value);
        } else {
          // Else display default protein
          this[key] = getProtein('nsp12 - RdRp');
          this.urlParams.set(key, 'nsp12 - RdRp');
        }
      } else if (key === 'ageRange') {
        // AgeRange is not being used currently so ignore
        return;
      } else if (key === 'customCoordinates' || key === 'residueCoordinates') {
        // If coordinates are specified, save them as an array of numbers
        // Coordinates can stay as a string in the URL
        let arr = [];
        value = encodeURIComponent(value);
        value = value.split('%2C');
        value.forEach((item, i) => {
          value[i] = parseInt(item);

          if (i % 2 === 1) {
            arr.push([value[i-1], value[i]]);
          }
        });

        this[key] = arr;
      } else if (key === 'tab') {
        // Check if the specified tab value is valid (included in TABS)
        if (Object.values(TABS).includes(value)) {
          this[key] = value;
        } else {
          // If not valid, set to home
          this.urlParams.set(key, TABS.TAB_EXAMPLE);
          this[key] = TABS.TAB_EXAMPLE;
        }
      } else if (key.includes('valid')) {
        // Convert the boolean config values to booleans
        this[key] = Boolean(value);
      } else {
        this[key] = value;
      }
    });

    // Update URL
    updateURLFromParams(this.urlParams);

    if (this.selectedLocationNodes.length == 0) {
      // If no locations in url, set default selected locations
      this.selectedLocationNodes = [
        getLocationByNameAndLevel(
          this.locationDataStoreInstance.selectTree,
          'USA',
          'country',
          true
        )[0],
        getLocationByNameAndLevel(
          this.locationDataStoreInstance.selectTree,
          'Canada',
          'country',
          true
        )[0],
      ].filter((node) => node !== undefined);
      initialConfigValues['selectedLocationNodes'] = this.selectedLocationNodes;
      this.initialConfigValues[
        'selectedLocationNodes'
      ] = this.selectedLocationNodes;
    }
  }

  // modifyQueryParams = autorun(() => {
  //   PARAMS_TO_TRACK.forEach((param) => {
  //     updateQueryStringParam(param, JSON.stringify(this[param]));
  //   });
  // });

  @action
  resetValues(values) {
    Object.keys(initialConfigValues).forEach((key) => {
      if (key in values) {
        this[key] = values[key];
      } else {
        this[key] = initialConfigValues[key];
      }

      // Special actions for some keys
      if (key === 'selectedLocationNodes') {
        this.locationDataStoreInstance.setSelectedNodes(values[key]);
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
    this.dataStoreInstance.fetchData();
  }

  @action
  applyPendingChanges(pending) {
    // Change table coloring settings when switching from DNA <-> AA
    if (this.dnaOrAa !== pending.dnaOrAa && pending.dnaOrAa === DNA_OR_AA.AA) {
      this.plotSettingsStoreInstance.tableColorMode =
        COLOR_MODES.COLOR_MODE_COMPARE;
      this.plotSettingsStoreInstance.tableCompareMode =
        COMPARE_MODES.COMPARE_MODE_MISMATCH;
      this.plotSettingsStoreInstance.tableCompareColor =
        COMPARE_COLORS.COLOR_MODE_ZAPPO;
    } else {
      // Clear table coloring settings
      this.plotSettingsStoreInstance.tableColorMode =
        COLOR_MODES.COLOR_MODE_COMPARE;
      this.plotSettingsStoreInstance.tableCompareMode =
        COMPARE_MODES.COMPARE_MODE_MISMATCH;
      this.plotSettingsStoreInstance.tableCompareColor =
        COMPARE_COLORS.COMPARE_COLOR_YELLOW;
    }

    // Overwrite any of our fields here with the pending ones
    Object.keys(pending).forEach((field) => {
      this[field] = pending[field];

      // Update urlParams
      this.urlParams.delete(field);

      if (field === 'selectedGene' || field === 'selectedProtein') {
        // Handle fields that return objects
        this.urlParams.set(field, pending[field].name);
      } else if (field === 'selectedMetadataFields') {
        // Ignore Metadata fields
        return;
      } else {
        this.urlParams.append(field, String(pending[field]));
      }
    });

    // selectedLocationNodes is displayed as node.level=node.value in the URL
    this.urlParams.delete('selectedLocationNodes');
    // selectedMetadataFields is displayed as key=value in the URL
    this.urlParams.delete('selectedMetadataFields');
    // ageRange is not currently being used
    this.urlParams.delete('ageRange');

    // Update the location node tree with our new selection
    this.locationDataStoreInstance.setSelectedNodes(this.selectedLocationNodes);

    // Update the URL params
    this.locationArray.forEach(name => {
      this.urlParams.delete(name);
    })

    this.selectedLocationNodes.forEach(node => {
      this.urlParams.append(String(node.level), String(node.value));
    })

    updateURLFromParams(this.urlParams);

    // Get the new data from the server
    this.dataStoreInstance.fetchData();
  }

  // Get a pretty name for the group
  getGroupLabel() {
    if (Object.keys(config.group_cols).includes(this.groupKey)) {
      return config.group_cols[this.groupKey].title;
    } else if (this.groupKey === GROUP_SNV) {
      if (this.dnaOrAa === DNA_OR_AA.DNA) {
        return 'NT SNV';
      } else {
        return 'AA SNV';
      }
    }
  }

  getCoordinateRanges() {
    // Set the coordinate range based off the coordinate mode
    if (this.coordinateMode === COORDINATE_MODES.COORD_GENE) {
      // Disable residue indices for non-protein-coding genes
      if (this.selectedGene.protein_coding === 0) {
        return this.selectedGene.ranges;
      }
      const coordinateRanges = [];
      this.residueCoordinates.forEach((range) => {
        // Make a deep copy of the current range
        const curRange = range.slice();
        for (let i = 0; i < this.selectedGene.aa_ranges.length; i++) {
          const curAARange = this.selectedGene.aa_ranges[i];
          const curNTRange = this.selectedGene.ranges[i];
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
          const curNTRange = this.selectedProtein.ranges[i];
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
  updateHoverGroup(group) {
    // console.log('UPDATE HOVER GROUP', group);
    if (group === null) {
      this.hoverGroup = null;
    } else {
      this.hoverGroup = group;
    }
  }

  @action
  updateSelectedGroups(groups) {
    this.selectedGroups = groups;
    if (this.groupKey === GROUP_SNV) {
      this.dataStoreInstance.processSelectedSnvs();
    }
  }

  getSelectedGroupIds() {
    const {
      dnaSnvMap,
      geneAaSnvMap,
      proteinAaSnvMap,
    } = this.snpDataStoreInstance;

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
    } = this.snpDataStoreInstance;

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
    } = this.snpDataStoreInstance;

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
  updateHoverLocation(location) {
    this.hoverLocation = location;
  }

  @action
  updateFocusedLocations(locations) {
    this.focusedLocations = locations;
  }
}

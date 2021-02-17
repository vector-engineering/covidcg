import {
  observable,
  action,
  toJS,
  //intercept, autorun
} from 'mobx';

import { getGene, getProtein } from '../utils/gene_protein';
import { queryReferenceSequence } from '../utils/reference';
import { getLocationByNameAndLevel } from '../utils/location';
import { intToISO, ISOToInt } from '../utils/date';

import {
  GROUP_SNV,
  DNA_OR_AA,
  COORDINATE_MODES,
  LOW_FREQ_FILTER_TYPES,
  MIN_DATE,
  COLOR_MODES,
  COMPARE_MODES,
  COMPARE_COLORS,
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

  dateRange: [-1, -1], // No initial date range
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
};

const urlParams = new URLSearchParams(window.location.search);

const defaultsFromParams = {};

PARAMS_TO_TRACK.forEach((param) => {
  // console.log('getting: ', param, urlParams.get(param));
  defaultsFromParams[param] = urlParams.get(param);
});

export class ConfigStore {
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
  @observable validCustomCoordinates = true;

  @observable customSequences = initialConfigValues.customSequences;
  @observable validCustomSequences = true;

  @observable residueCoordinates = initialConfigValues.residueCoordinates;
  @observable validResidueCoordinates = true;

  @observable coordinateMode = initialConfigValues.coordinateMode;

  @observable dateRange = initialConfigValues.dateRange;
  @observable startDate = initialConfigValues.startDate;
  @observable endDate = initialConfigValues.endDate;
  @observable validDateRange = true;

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

  constructor() {}

  init() {
    this.plotSettingsStoreInstance = rootStoreInstance.plotSettingsStore;
    this.locationDataStoreInstance = rootStoreInstance.locationDataStore;
    this.dataStoreInstance = rootStoreInstance.dataStore;

    PARAMS_TO_TRACK.forEach((param) => {
      if (defaultsFromParams[param]) {
        // console.log('setting: ', param, urlParams.get(param));
        this[param] = defaultsFromParams[param];
      }
    });

    // Set default selected locations
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
    this.dataStoreInstance.updateCaseData(() => {});
  }

  @action
  changeGrouping(groupKey, dnaOrAa) {
    if (this.groupKey !== groupKey) {
      // If groupings were changed, then clear selected groups
      this.selectedGroups = [];
    } else if (groupKey === GROUP_SNV && this.dnaOrAa !== dnaOrAa) {
      // While in SNV mode, if we switched from DNA to AA, or vice-versa,
      // then clear selected groups
      this.selectedGroups = [];
    }

    // Change table coloring settings when switching from DNA <-> AA
    if (this.dnaOrAa !== dnaOrAa && dnaOrAa === DNA_OR_AA.AA) {
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

    this.groupKey = groupKey;
    this.dnaOrAa = dnaOrAa;

    // If we switched to non-SNP grouping in AA-mode,
    // then make sure we don't have "All Genes" or "All Proteins" selected
    if (this.groupKey !== GROUP_SNV && this.dnaOrAa === DNA_OR_AA.AA) {
      if (this.selectedGene.name === 'All Genes') {
        // Switch back to S gene
        this.selectedGene = getGene('S');
      }
      if (this.selectedProtein.name === 'All Proteins') {
        // Switch back to nsp12 protein
        this.selectedProtein = getProtein('nsp12 - RdRp');
      }
    }

    // this.dataStoreInstance.updateCaseData();
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

  setDefaultGeneResidueCoordinates() {
    if (this.selectedGene.name === 'All Genes') {
      this.residueCoordinates = [];
    } else {
      this.residueCoordinates = [[1, this.selectedGene.len_aa]];
    }
    // Clear valid residue coordinates flag
    this.validResidueCoordinates = true;
  }

  setDefaultProteinResidueCoordinates() {
    if (this.selectedProtein.name === 'All Proteins') {
      this.residueCoordinates = [];
    } else {
      this.residueCoordinates = [[1, this.selectedProtein.len_aa]];
    }
    // Clear valid residue coordinates flag
    this.validResidueCoordinates = true;
  }

  @action
  updateCoordinateMode(coordinateMode) {
    this.coordinateMode = coordinateMode;

    // If we switched to a coordinate mode that doesn't support AA SNPs,
    // then switch off of it now
    if (
      this.dnaOrAa === DNA_OR_AA.AA &&
      this.coordinateMode !== COORDINATE_MODES.COORD_GENE &&
      this.coordinateMode !== COORDINATE_MODES.COORD_PROTEIN
    ) {
      this.dnaOrAa = DNA_OR_AA.DNA;
    }

    if (this.coordinateMode === COORDINATE_MODES.COORD_GENE) {
      this.setDefaultGeneResidueCoordinates();
    } else if (this.coordinateMode === COORDINATE_MODES.COORD_PROTEIN) {
      this.setDefaultProteinResidueCoordinates();
    }
  }

  @action
  updateSelectedGene(selectedGene) {
    const initialSelectedGene = this.selectedGene;
    this.selectedGene = getGene(selectedGene);

    // If we selected a new gene, then update the residue coordinates
    if (this.selectedGene.name !== initialSelectedGene.name) {
      this.setDefaultGeneResidueCoordinates();
    }
  }

  @action
  updateSelectedProtein(selectedProtein) {
    const initialSelectedProtein = this.selectedProtein;
    this.selectedProtein = getProtein(selectedProtein);

    // If we selected a new protein, then update the residue coordinates
    if (this.selectedProtein.name !== initialSelectedProtein.name) {
      this.setDefaultProteinResidueCoordinates();
    }
  }

  @action
  updateResidueCoordinates(residueCoordinates) {
    this.residueCoordinates = residueCoordinates;
  }
  @action
  updateValidResidueCoordinates(valid) {
    this.validResidueCoordinates = valid;
  }

  @action
  updateSelectedPrimers(selectedPrimers) {
    this.selectedPrimers = selectedPrimers;
  }

  @action
  updateCustomCoordinates(customCoordinates) {
    this.customCoordinates = customCoordinates;
  }
  @action
  updateValidCustomCoordinates(valid) {
    this.validCustomCoordinates = valid;
  }

  @action
  updateCustomSequences(customSequences) {
    this.customSequences = customSequences;
  }
  @action
  updateValidCustomSequences(valid) {
    this.validCustomSequences = valid;
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

  @action
  selectLocations(selectedLocationNodes) {
    this.selectedLocationNodes = selectedLocationNodes;

    // Clear metadata fields
    this.selectedMetadataFields = {};

    // if (!selectedLocationNodes || !selectedLocationNodes[0]) {
    //   this.locationDataStoreInstance.deselectAll();
    //   this.dataStoreInstance.emptyCaseData();
    // } else {
    //   this.dataStoreInstance.updateCaseData();
    // }
  }

  @action
  updateSelectedMetadataFields(field, options) {
    this.selectedMetadataFields[field] = options;
    // this.dataStoreInstance.updateCaseData();
  }

  @action
  updateAgeRange(ageRange) {
    this.ageRange = ageRange;
  }

  @action
  updateDateRange(dateRange) {
    this.dateRange = dateRange;
    // this.dataStoreInstance.updateAggCaseDataByGroup();
    // if (this.groupKey === GROUP_SNV) {
    //   this.dataStoreInstance.processCooccurrenceData();
    // }
  }
  @action
  updateStartDate(startDate) {
    this.startDate = startDate;
    this.checkValidDateRanges();
  }
  @action
  updateEndDate(endDate) {
    this.endDate = endDate;
    this.checkValidDateRanges();
  }
  @action
  updateValidDateRange(validDateRange) {
    this.validDateRange = validDateRange;
  }
  checkValidDateRanges() {
    const validDateRange = ISOToInt(this.startDate) < ISOToInt(this.endDate);
    this.updateValidDateRange(validDateRange);
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
    this.dataStoreInstance.processSelectedSnvs();
  }

  @action
  updateHoverLocation(location) {
    this.hoverLocation = location;
  }

  @action
  updateFocusedLocations(locations) {
    this.focusedLocations = locations;
  }

  @action
  setLowFreqFilterType(lowFreqFilterType) {
    this.lowFreqFilterType = lowFreqFilterType;
  }

  @action
  setMinLocalCounts(minLocalCounts) {
    this.minLocalCounts = minLocalCounts;
  }

  @action
  setMinGlobalCounts(minGlobalCounts) {
    this.minGlobalCounts = minGlobalCounts;
  }

  @action
  setMaxGroupCounts(maxGroupCounts) {
    this.maxGroupCounts = maxGroupCounts;
  }
}

import {
  observable,
  action,
  toJS,
  //intercept, autorun
} from 'mobx';
import { dataStoreInstance, plotSettingsStoreInstance } from './rootStore';
import _ from 'underscore';

import { getGene, getProtein } from '../utils/gene_protein';
import {
  getLocationByNameAndLevel,
  loadSelectTree,
  assignObjectPaths,
  getNodeFromPath,
  deselectAll,
  selectAll,
} from '../utils/location';
import { queryReferenceSequence } from '../utils/reference';

import {
  COLOR_MODES,
  COMPARE_MODES,
  COMPARE_COLORS,
} from '../constants/plotSettings';
import {
  GROUP_KEYS,
  DNA_OR_AA,
  COORDINATE_MODES,
  LOW_FREQ_FILTER_TYPES,
} from '../constants/config';

// import { updateQueryStringParam } from '../utils/updateQueryParam';
import { PARAMS_TO_TRACK } from './paramsToTrack';

// Define initial values

const initialSelectTree = Object.assign(loadSelectTree(), {
  expanded: true,
});
assignObjectPaths(initialSelectTree);

export const initialConfigValues = {
  groupKey: GROUP_KEYS.GROUP_SNV,
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

  selectTree: initialSelectTree,
  selectedLocationNodes: [
    getLocationByNameAndLevel(initialSelectTree, 'USA', 'country', true)[0],
    getLocationByNameAndLevel(initialSelectTree, 'Canada', 'country', true)[0],
  ],

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
  minLocalCountsToShow: 50,
  minGlobalCountsToShow: 100,
};

const urlParams = new URLSearchParams(window.location.search);

const defaultsFromParams = {};

PARAMS_TO_TRACK.forEach((param) => {
  // console.log('getting: ', param, urlParams.get(param));
  defaultsFromParams[param] = urlParams.get(param);
});

class ObservableConfigStore {
  @observable groupKey = initialConfigValues.groupKey;
  @observable dnaOrAa = initialConfigValues.dnaOrAa;

  @observable selectedGene = initialConfigValues.selectedGene;
  @observable selectedProtein = initialConfigValues.selectedProtein;
  @observable selectedPrimers = initialConfigValues.selectedPrimers;
  @observable customCoordinates = initialConfigValues.customCoordinates;
  @observable customSequences = initialConfigValues.customSequences;
  @observable residueCoordinates = initialConfigValues.residueCoordinates;

  @observable coordinateMode = initialConfigValues.coordinateMode;

  @observable dateRange = initialConfigValues.dateRange;

  @observable selectTree = initialConfigValues.selectTree;
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
  @observable minLocalCountsToShow = initialConfigValues.minLocalCountsToShow;
  @observable minGlobalCountsToShow = initialConfigValues.minGlobalCountsToShow;

  constructor() {
    PARAMS_TO_TRACK.forEach((param) => {
      if (defaultsFromParams[param]) {
        // console.log('setting: ', param, urlParams.get(param));
        this[param] = defaultsFromParams[param];
      }
    });
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
        // Get the current tree
        const selectTree = Object.assign({}, this.selectTree);
        deselectAll(selectTree);

        // Select the specified nodes
        values[key].forEach((node) => {
          const nodeObj =
            'path' in node
              ? getNodeFromPath(selectTree, node['path'])
              : selectTree;
          // Select all of the nodes children
          selectAll(nodeObj);
        });

        this.selectTree = selectTree;
      }
    });

    // Manually set residue coordinates if they weren't specified
    if (
      'selectedGene' in values &&
      !('residueCoordinates' in values) &&
      this.selectedGene.gene !== 'All Genes'
    ) {
      this.residueCoordinates = [[1, this.selectedGene.len_aa]];
    } else if (
      'selectedProtein' in values &&
      !('residueCoordinates' in values) &&
      this.selectedProtein.protein !== 'All Proteins'
    ) {
      this.residueCoordinates = [[1, this.selectedProtein.len_aa]];
    }

    // Trigger data re-run
    dataStoreInstance.updateCaseData(() => {});
  }

  @action
  changeGrouping(groupKey, dnaOrAa) {
    // console.log(
    //   'CHANGE GROUPING. GROUP KEY:',
    //   groupKey,
    //   'DNA OR AA:',
    //   dnaOrAa
    // );

    if (this.groupKey !== groupKey) {
      // If groupings were changed, then clear selected groups
      this.selectedGroups = [];
    } else if (groupKey === GROUP_KEYS.GROUP_SNV && this.dnaOrAa !== dnaOrAa) {
      // While in SNV mode, if we switched from DNA to AA, or vice-versa,
      // then clear selected groups
      this.selectedGroups = [];
    }

    // Change table coloring settings when switching from DNA <-> AA
    if (this.dnaOrAa !== dnaOrAa && dnaOrAa === DNA_OR_AA.AA) {
      plotSettingsStoreInstance.tableColorMode = COLOR_MODES.COLOR_MODE_COMPARE;
      plotSettingsStoreInstance.tableCompareMode =
        COMPARE_MODES.COMPARE_MODE_MISMATCH;
      plotSettingsStoreInstance.tableCompareColor =
        COMPARE_COLORS.COLOR_MODE_ZAPPO;
    } else {
      // Clear table coloring settings
      plotSettingsStoreInstance.tableColorMode = COLOR_MODES.COLOR_MODE_COMPARE;
      plotSettingsStoreInstance.tableCompareMode =
        COMPARE_MODES.COMPARE_MODE_MISMATCH;
      plotSettingsStoreInstance.tableCompareColor =
        COMPARE_COLORS.COMPARE_COLOR_YELLOW;
    }

    this.groupKey = groupKey;
    this.dnaOrAa = dnaOrAa;

    // If we switched to non-SNP grouping in AA-mode,
    // then make sure we don't have "All Genes" or "All Proteins" selected
    if (
      this.groupKey !== GROUP_KEYS.GROUP_SNV &&
      this.dnaOrAa === DNA_OR_AA.AA
    ) {
      if (this.selectedGene.gene === 'All Genes') {
        // Switch back to S gene
        this.selectedGene = getGene('S');
      }
      if (this.selectedProtein.protein === 'All Proteins') {
        // Switch back to nsp12 protein
        this.selectedProtein = getProtein('nsp12 - RdRp');
      }
    }

    dataStoreInstance.updateCaseData();
  }

  // Get a pretty name for the group
  getGroupLabel() {
    if (this.groupKey === GROUP_KEYS.GROUP_LINEAGE) {
      return 'Lineage';
    } else if (this.groupKey === GROUP_KEYS.GROUP_CLADE) {
      return 'Clade';
    } else if (this.groupKey === GROUP_KEYS.GROUP_SNV) {
      if (this.dnaOrAa === DNA_OR_AA.DNA) {
        return 'NT SNV';
      } else {
        return 'AA SNV';
      }
    }
  }

  @action
  changeCoordinateMode({
    coordinateMode,
    selectedGene,
    selectedProtein,
    residueCoordinates,
    selectedPrimers,
    customCoordinates,
    customSequences,
  }) {
    let initial = Object.assign({
      coordinateMode: toJS(this.coordinateMode),
      selectedGene: toJS(this.selectedGene),
      selectedProtein: toJS(this.selectedProtein),
      residueCoordinates: toJS(this.residueCoordinates),
      selectedPrimers: toJS(this.selectedPrimers),
      customCoordinates: toJS(this.customCoordinates),
      customSequences: toJS(this.customSequences),
    });
    // console.log(initial);

    this.coordinateMode = coordinateMode;
    this.selectedGene = getGene(selectedGene);
    this.selectedProtein = getProtein(selectedProtein);
    this.residueCoordinates = residueCoordinates;
    this.selectedPrimers = selectedPrimers;
    this.customCoordinates = customCoordinates;
    this.customSequences = customSequences;

    // If we selected a new gene/protein, then update the residue coordinates
    if (
      this.coordinateMode === COORDINATE_MODES.COORD_GENE &&
      (this.selectedGene.gene !== initial.selectedGene.gene ||
        this.coordinateMode !== initial.coordinateMode)
    ) {
      if (this.selectedGene.gene === 'All Genes') {
        this.residueCoordinates = [];
      } else {
        this.residueCoordinates = [[1, this.selectedGene.len_aa]];
      }
    } else if (
      this.coordinateMode === COORDINATE_MODES.COORD_PROTEIN &&
      (this.selectedProtein.protein !== initial.selectedProtein.protein ||
        this.coordinateMode !== initial.coordinateMode)
    ) {
      if (this.selectedProtein.protein === 'All Proteins') {
        this.residueCoordinates = [];
      } else {
        this.residueCoordinates = [[1, this.selectedProtein.len_aa]];
      }
    }

    // If we switched to a coordinate mode that doesn't support AA SNPs,
    // then switch off of it now
    if (
      this.dnaOrAa === DNA_OR_AA.AA &&
      this.coordinateMode !== COORDINATE_MODES.COORD_GENE &&
      this.coordinateMode !== COORDINATE_MODES.COORD_PROTEIN
    ) {
      this.dnaOrAa = DNA_OR_AA.DNA;
    }

    // If nothing changed, then skip the update
    if (this.coordinateMode !== initial.coordinateMode) {
      // Do nothing
    } else if (
      this.coordinateMode === COORDINATE_MODES.COORD_GENE &&
      this.selectedGene.gene === initial.selectedGene.gene &&
      _.isEqual(toJS(this.residueCoordinates), initial.residueCoordinates)
    ) {
      return;
    } else if (
      this.coordinateMode === COORDINATE_MODES.COORD_PROTEIN &&
      this.selectedProtein.protein == initial.selectedProtein.protein &&
      _.isEqual(toJS(this.residueCoordinates), initial.residueCoordinates)
    ) {
      return;
    } else if (this.coordinateMode === COORDINATE_MODES.COORD_PRIMER) {
      let changed = false;
      if (this.selectedPrimers.length !== initial.selectedPrimers.length) {
        changed = true;
      } else {
        for (let i = 0; i < this.selectedPrimers.length; i++) {
          if (!_.isEqual(this.selectedPrimers[i], initial.selectedPrimers[i])) {
            changed = true;
            break;
          }
        }
      }
      if (!changed) {
        return;
      }
    } else if (
      this.coordinateMode === COORDINATE_MODES.COORD_CUSTOM &&
      _.isEqual(toJS(this.customCoordinates), initial.customCoordinates)
    ) {
      return;
    } else if (
      this.coordinateMode === COORDINATE_MODES.COORD_SEQUENCE &&
      _.isEqual(toJS(this.customSequences), initial.customSequences)
    ) {
      return;
    }

    // Clear selected groups?
    // TODO: we don't need to do this, depending on the selection.
    //       do this in a smarter way
    this.selectedGroups = [];

    dataStoreInstance.updateCaseData();
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
          if (curRange[0] >= curAARange[0] && curRange[0] <= curAARange[1]) {
            coordinateRanges.push([
              curNTRange[0] + (curRange[0] - 1) * 3,
              curNTRange[0] - 1 + Math.min(curRange[1], curAARange[1]) * 3,
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
          if (curRange[0] >= curAARange[0] && curRange[0] <= curAARange[1]) {
            coordinateRanges.push([
              curNTRange[0] + (curRange[0] - 1) * 3,
              curNTRange[0] - 1 + Math.min(curRange[1], curAARange[1]) * 3,
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

    if (!selectedLocationNodes || !selectedLocationNodes[0]) {
      this.selectTree = deselectAll(toJS(this.selectTree));
      dataStoreInstance.emptyCaseData();
    } else {
      dataStoreInstance.updateCaseData();
    }
  }

  @action
  updateSelectedMetadataFields(selectedMetadataFields, ageRange) {
    this.selectedMetadataFields = selectedMetadataFields;
    this.ageRange = ageRange;
    dataStoreInstance.updateCaseData();
  }

  @action
  selectDateRange(dateRange) {
    this.dateRange = dateRange;
    dataStoreInstance.updateAggCaseDataByGroup();
    if (this.groupKey === GROUP_KEYS.GROUP_SNV) {
      dataStoreInstance.processCooccurrenceData();
    }
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
    dataStoreInstance.processSelectedSnvs();
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
  setLowFreqFilters({
    lowFreqFilterType,
    minLocalCountsToShow,
    minGlobalCountsToShow,
    maxGroupCounts,
  }) {
    this.lowFreqFilterType = lowFreqFilterType;
    this.minLocalCountsToShow = minLocalCountsToShow;
    this.minGlobalCountsToShow = minGlobalCountsToShow;
    this.maxGroupCounts = maxGroupCounts;
    dataStoreInstance.updateCaseData();
  }
}

export default ObservableConfigStore;

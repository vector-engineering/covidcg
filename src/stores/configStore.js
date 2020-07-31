import { observable, action, toJS } from 'mobx';
import { dataStoreInstance } from './rootStore';
import _ from 'underscore';

import { getGene } from '../utils/gene';
import { getProtein } from '../utils/protein';
import {
  getLocationIds,
  getLocationByNameAndLevel,
  loadSelectTree,
} from '../utils/location';

const LOCAL_COUNTS = 'LOCAL_COUNTS';
const GLOBAL_COUNTS = 'GLOBAL_COUNTS';
const GROUP_COUNTS = 'GROUP_COUNTS';

export const LOW_FREQ_FILTER_TYPES = {
  LOCAL_COUNTS,
  GLOBAL_COUNTS,
  GROUP_COUNTS,
};

class ObservableConfigStore {
  @observable groupKey = 'lineage'; // lineage, clade, snp
  @observable dnaOrAa = 'dna';

  // COORDINATE SETTINGS
  // Select the Spike gene and nsp13 protein by default
  @observable selectedGene = getGene('S');
  @observable selectedProtein = getProtein('nsp13');
  @observable selectedPrimers = [];
  @observable customCoordinates = [8000, 12000];

  // Selecting the gene as the coordinate range by default
  @observable coordinateMode = 'gene';
  @observable coordinateRanges = this.selectedGene.ranges;

  @observable dateRange = [-1, -1]; // No initial date range

  @observable selectedLocationNodes = [];
  @observable selectedLocationIds = [];

  @observable hoverGroup = null;
  // @observable selectedGroups = [{ group: 'B.1' }, { group: 'B.1.3' }];
  @observable selectedGroups = [];

  // Metadata filtering
  @observable selectedMetadataFields = {};
  @observable ageRange = [null, null];

  // Location tab
  @observable hoverLocation = null;
  @observable focusedLocations = []; // Selected locations in the location tab

  @observable lowFreqFilterType = LOCAL_COUNTS;
  @observable maxLineagesToShow = 10;
  @observable minLocalCountsToShow = 50;

  constructor() {
    const selectTree = loadSelectTree();
    // Select NYC by default
    let NYCNode = getLocationByNameAndLevel(
      selectTree,
      'New York City',
      'location'
    );
    NYCNode[0].checked = true;
    let NYCLocationId = getLocationIds(NYCNode);

    let MassNode = getLocationByNameAndLevel(
      selectTree,
      'Massachusetts',
      'division'
    );
    MassNode[0].checked = true;
    let MassLocationId = getLocationIds(MassNode);

    this.selectedLocationIds = NYCLocationId.concat(MassLocationId);
    // console.log(this.selectedLocationIds);
    this.selectedLocationNodes = [NYCNode[0], MassNode[0]];
  }

  @action
  changeGrouping(groupKey, dnaOrAa) {
    // console.log(
    //   'CHANGE GROUPING. GROUP KEY:',
    //   groupKey,
    //   'DNA OR AA:',
    //   dnaOrAa
    // );

    // Clear selected groups?
    if (this.groupKey !== groupKey) {
      this.selectedGroups = [];
    } else if (groupKey === 'snp' && this.dnaOrAa !== dnaOrAa) {
      this.selectedGroups = [];
    }

    this.groupKey = groupKey;
    this.dnaOrAa = dnaOrAa;

    // If we switched to non-SNP grouping in AA-mode,
    // then make sure we don't have "All Genes" or "All Proteins" selected
    if (this.groupKey !== 'snp' && this.dnaOrAa === 'aa') {
      if (this.selectedGene.gene === 'All Genes') {
        // Switch back to S gene
        this.selectedGene = getGene('S');
      }
      if (this.selectedProtein.protein === 'All Proteins') {
        // Switch back to nsp13 protein
        this.selectedProtein = getProtein('nsp13');
      }
    }

    dataStoreInstance.updateCaseData();
  }

  // Get a pretty name for the group
  getGroupLabel() {
    if (this.groupKey === 'lineage') {
      return 'Lineage';
    } else if (this.groupKey === 'clade') {
      return 'Clade';
    } else if (this.groupKey === 'snp') {
      if (this.dnaOrAa === 'dna') {
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
    selectedPrimers,
    customCoordinates,
  }) {
    // console.log('CHANGE COORDINATE MODE', coordinateMode);
    // console.log('SELECTED GENE:', selectedGene);
    // console.log('SELECTED PROTEIN:', selectedProtein);
    // console.log('SELECTED PRIMERS:', selectedPrimers);
    // console.log('CUSTOM COORDINATES:', customCoordinates);

    let initial = Object.assign({
      coordinateMode: toJS(this.coordinateMode),
      selectedGene: toJS(this.selectedGene),
      selectedProtein: toJS(this.selectedProtein),
      selectedPrimers: toJS(this.selectedPrimers),
      customCoordinates: toJS(this.customCoordinates),
    });
    // console.log(initial);

    this.coordinateMode = coordinateMode;
    this.selectedGene = getGene(selectedGene);
    this.selectedProtein = getProtein(selectedProtein);
    this.selectedPrimers = selectedPrimers;
    this.customCoordinates = customCoordinates;

    // Set the coordinate range based off the coordinate mode
    if (coordinateMode === 'gene') {
      this.coordinateRanges = this.selectedGene.ranges;
    } else if (coordinateMode === 'protein') {
      this.coordinateRanges = this.selectedProtein.ranges;
    } else if (coordinateMode === 'primer') {
      let ranges = [];
      this.selectedPrimers.forEach((primer) => {
        ranges.push([primer.Start, primer.End]);
      });
      this.coordinateRanges = ranges;
    } else if (coordinateMode === 'custom') {
      this.coordinateRanges = [this.customCoordinates];
    }

    // If we switched to a coordinate mode that doesn't support AA SNPs,
    // then switch off of it now
    if (
      this.dnaOrAa === 'aa' &&
      this.coordinateMode !== 'gene' &&
      this.coordinateMode !== 'protein'
    ) {
      this.dnaOrAa = 'dna';
    }

    // If nothing changed, then skip the update
    if (this.coordinateMode !== initial.coordinateMode) {
      // Do nothing
    } else if (
      this.coordinateMode === 'gene' &&
      this.selectedGene.gene === initial.selectedGene.gene
    ) {
      return;
    } else if (
      this.coordinateMode === 'protein' &&
      this.selectedProtein.protein == initial.selectedProtein.protein
    ) {
      return;
    } else if (this.coordinateMode === 'primer') {
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
      this.coordinateMode === 'custom' &&
      this.coordinateRanges[0][0] === initial.customCoordinates[0][0] &&
      this.coordinateRanges[0][1] === initial.customCoordinates[0][1]
    ) {
      return;
    }

    // Clear selected groups?
    // TODO: we don't need to do this, depending on the selection.
    //       do this in a smarter way
    this.selectedGroups = [];

    dataStoreInstance.updateCaseData();
  }

  @action
  selectLocations(selectedNodes) {
    this.selectedLocationNodes = selectedNodes;
    this.selectedLocationIds = getLocationIds(selectedNodes);

    // Clear metadata fields
    this.selectedMetadataFields = {};

    if (!selectedNodes || !selectedNodes[0]) {
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
  }

  @action
  updateHoverGroup(group) {
    // console.log('UPDATE HOVER GROUP', group);
    if (group === null) {
      this.hoverGroup = null;
    } else if (!dataStoreInstance.groupsToKeep.includes(group)) {
      this.hoverGroup = 'other';
    } else {
      this.hoverGroup = group;
    }
  }

  @action
  updateSelectedGroups(groups) {
    this.selectedGroups = groups;
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
  setLowFreqFilterType(type) {
    this.lowFreqFilterType = type;
    dataStoreInstance.updateGroupsToKeep();
  }

  @action
  setMaxLineages(num) {
    this.maxLineagesToShow = num;
    dataStoreInstance.updateGroupsToKeep();
  }

  @action
  setMinLocalCounts(num) {
    this.minLocalCountsToShow = num;
    dataStoreInstance.updateGroupsToKeep();
  }
}

export default ObservableConfigStore;

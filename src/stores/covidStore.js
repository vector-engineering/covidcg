import { observable, action, toJS } from 'mobx';
import _ from 'underscore';

import {
  processCaseData,
  aggCaseDataByGroup,
} from '../utils/caseDataWorkerWrapper';
import {
  downloadAcknowledgements,
  downloadAggCaseData,
} from '../utils/downloadWorkerWrapper';
import { getGene } from '../utils/gene';
import { getProtein } from '../utils/protein';
//import { getLineagesFromGene } from '../utils/lineageData';
import {
  loadSelectTree,
  getLocationByNameAndLevel,
  getLocationIds,
} from '../utils/location';
import { downloadBlobURL, generateSelectionString } from '../utils/download';
import { uiStoreInstance } from './rootStore';

class ObservableCovidStore {
  @observable groupKey = 'lineage';
  @observable dnaOrAa = 'dna';

  // Select the Spike gene and nsp13 protein by default
  @observable selectedGene = getGene('S');
  @observable selectedProtein = getProtein('nsp13');
  @observable selectedPrimers = [];
  @observable customCoordinates = [8000, 12000];

  // Selecting the gene as the coordinate range by default
  @observable coordinateMode = 'gene';
  @observable coordinateRanges = [
    [this.selectedGene.start, this.selectedGene.end],
  ];

  @observable selectTree = loadSelectTree();
  @observable selectedLocationIds = [];
  @observable caseData = [];
  @observable changingPositions = {};
  @observable caseDataAggGroup = [];
  @observable dateRange = [-1, -1]; // No initial date range
  @observable selectedRows = [];
  @observable groupsToKeep = {};

  @observable hoverGroup = null;
  @observable selectedGroups = [];

  // Metadata filtering
  @observable numSequencesBeforeMetadataFiltering = 0;
  @observable metadataCounts = {};
  @observable selectedMetadataFields = {};
  @observable ageRange = [null, null];

  constructor() {
    // Select NYC by default
    let NYCNode = getLocationByNameAndLevel(
      this.selectTree,
      'New York City',
      'location'
    );
    NYCNode[0].checked = true;
    let NYCLocationId = getLocationIds(NYCNode);
    this.selectedLocationIds = NYCLocationId;

    uiStoreInstance.onCaseDataStateStarted();

    this.updateCaseData();
  }

  @action
  changeGrouping(_groupKey, _dnaOrAa) {
    console.log(
      'CHANGE GROUPING. GROUP KEY:',
      _groupKey,
      'DNA OR AA:',
      _dnaOrAa
    );
    this.groupKey = _groupKey;
    this.dnaOrAa = _dnaOrAa;

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

    this.updateCaseData();
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
      this.coordinateRanges = [
        [this.selectedGene.start, this.selectedGene.end],
      ];
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

    this.updateCaseData();
  }

  @action
  selectLocations(selectedNodes) {
    this.selectedLocationIds = getLocationIds(selectedNodes);

    if (!selectedNodes || !selectedNodes[0]) {
      this.emptyCaseData();
    } else {
      this.updateCaseData();
    }
  }

  @action
  updateSelectedMetadataFields(selectedMetadataFields, ageRange) {
    this.selectedMetadataFields = selectedMetadataFields;
    this.ageRange = ageRange;
    this.updateCaseData();
  }

  @action
  selectDateRange(_dateRange) {
    this.dateRange = _dateRange;
    this.updateAggCaseDataByGroup();
  }

  @action
  updateAggCaseDataByGroup(suppressUIUpdate = false) {
    suppressUIUpdate ? null : uiStoreInstance.onAggCaseDataStarted();
    aggCaseDataByGroup(
      {
        caseData: toJS(this.caseData),
        coordinateMode: toJS(this.coordinateMode),
        coordinateRanges: toJS(this.coordinateRanges),
        selectedGene: toJS(this.selectedGene),
        selectedProtein: toJS(this.selectedProtein),
        groupKey: toJS(this.groupKey),
        dnaOrAa: toJS(this.dnaOrAa),
        dateRange: toJS(this.dateRange),
      },
      ({ caseDataAggGroup, changingPositions, groupsToKeepObj }) => {
        // console.log(caseDataAggGroup);
        this.caseDataAggGroup = caseDataAggGroup;
        this.changingPositions = changingPositions;
        this.groupsToKeep = groupsToKeepObj;
        console.log('AGG_CASE_DATA FINISHED');

        suppressUIUpdate ? null : uiStoreInstance.onAggCaseDataFinished();
        suppressUIUpdate ? null : uiStoreInstance.onCaseDataStateFinished();
      }
    );
  }

  @action
  emptyCaseData() {
    this.caseData = [];
    this.selectedRows = [];
    this.caseDataAggGroup = [];
    this.changingPositions = {};
  }

  @action
  updateCaseData(suppressUIUpdate = false) {
    suppressUIUpdate ? null : uiStoreInstance.onCaseDataStateStarted();

    processCaseData(
      {
        selectedLocationIds: toJS(this.selectedLocationIds),
        coordinateMode: toJS(this.coordinateMode),
        coordinateRanges: toJS(this.coordinateRanges),
        selectedGene: toJS(this.selectedGene),
        selectedProtein: toJS(this.selectedProtein),
        groupKey: toJS(this.groupKey),
        dnaOrAa: toJS(this.dnaOrAa),
        selectedMetadataFields: toJS(this.selectedMetadataFields),
        ageRange: toJS(this.ageRange),
      },
      ({
        aggCaseDataList,
        selectedRows,
        metadataCounts,
        numSequencesBeforeMetadataFiltering,
      }) => {
        this.caseData = aggCaseDataList;
        this.selectedRows = selectedRows;
        this.metadataCounts = metadataCounts;
        this.numSequencesBeforeMetadataFiltering = numSequencesBeforeMetadataFiltering;
        console.log('CASE_DATA FINISHED');

        this.updateAggCaseDataByGroup((suppressUIUpdate = false));
      }
    );
  }

  @action
  updateHoverGroup(group) {
    // console.log('UPDATE HOVER GROUP', group);
    this.hoverGroup = group;
  }

  @action
  updateSelectedGroups(groups) {
    this.selectedGroups = groups;
  }

  @action
  downloadAcknowledgements() {
    // console.log('DOWNLOAD ACKNOWLEDGEMENTS');
    downloadAcknowledgements(
      { selectedRows: toJS(this.selectedRows) },
      (res) => {
        // console.log(res);
        downloadBlobURL(
          res.blobURL,
          generateSelectionString(
            'acknowledgements',
            'csv',
            this.groupKey,
            this.dnaOrAa,
            this.selectedGene,
            this.selectedLocationIds,
            this.dateRange
          )
        );
      }
    );
  }

  @action
  downloadAggCaseData() {
    downloadAggCaseData(
      {
        groupKey: this.groupKey,
        dnaOrAa: this.dnaOrAa,
        coordinateMode: this.coordinateMode,
        caseDataAggGroup: toJS(this.caseDataAggGroup),
      },
      (res) => {
        downloadBlobURL(
          res.blobURL,
          generateSelectionString(
            'agg_data',
            'csv',
            this.groupKey,
            this.dnaOrAa,
            this.selectedGene,
            this.selectedLocationIds,
            this.dateRange
          )
        );
      }
    );
  }
}

export default ObservableCovidStore;

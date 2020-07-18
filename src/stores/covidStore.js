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
  @observable groupKey = null;
  @observable dnaOrAa = null;

  @observable selectedGene = {};
  @observable selectedProtein = {};
  @observable selectedPrimers = [];
  @observable customCoordinates = [8000, 12000];

  @observable coordinateMode = null;
  @observable coordinateRanges = [];

  @observable selectTree = [];
  @observable selectedLocationIds = [];
  @observable caseData = [];
  @observable changingPositions = {};
  @observable caseDataAggGroup = [];
  @observable dateRange = [];
  @observable selectedRows = [];
  @observable groupsToKeep = {};

  constructor() {
    // Select the Spike gene and nsp13 protein by default
    let defaultGene = getGene('S');
    let defaultProtein = getProtein('nsp13');
    let defaultSelectedPrimers = [];
    let defaultCustomCoordinates = [8000, 12000];
    // Selecting the gene as the coordinate range by default
    let defaultCoordinateMode = 'gene';
    let defaultCoordinateRanges = [[defaultGene.start, defaultGene.end]];

    let selectTree = loadSelectTree();
    // Select NYC by default
    let NYCNode = getLocationByNameAndLevel(
      selectTree,
      'New York City',
      'location'
    );
    NYCNode[0].checked = true;
    let NYCLocationId = getLocationIds(NYCNode);

    let initialGroupKey = 'lineage';
    let initialDnaOrAa = 'dna';
    let initialLocationIds = NYCLocationId;

    // No initial date range
    let initialDateRange = [-1, -1];

    uiStoreInstance.onCaseDataStateStarted();

    processCaseData(
      {
        selectedLocationIds: toJS(initialLocationIds),
        coordinateMode: toJS(defaultCoordinateMode),
        coordinateRanges: toJS(defaultCoordinateRanges),
        selectedGene: toJS(defaultGene),
        selectedProtein: toJS(defaultProtein),
        groupKey: toJS(initialGroupKey),
        dnaOrAa: toJS(initialDnaOrAa),
      },
      ({ aggCaseDataList, selectedRows }) => {
        uiStoreInstance.onAggCaseDataStarted();

        aggCaseDataByGroup(
          {
            caseData: toJS(aggCaseDataList),
            coordinateMode: toJS(defaultCoordinateMode),
            coordinateRanges: toJS(defaultCoordinateRanges),
            selectedGene: toJS(defaultGene),
            selectedProtein: toJS(defaultProtein),
            groupKey: toJS(initialGroupKey),
            dnaOrAa: toJS(initialDnaOrAa),
            dateRange: toJS(initialDateRange),
          },
          ({ changingPositions, caseDataAggGroup, groupsToKeepObj }) => {
            this.groupsToKeep = groupsToKeepObj;
            this.groupKey = initialGroupKey;
            this.dnaOrAa = initialDnaOrAa;

            this.selectedGene = defaultGene;
            this.selectedProtein = defaultProtein;
            this.selectedPrimers = defaultSelectedPrimers;
            this.customCoordinates = defaultCustomCoordinates;

            this.coordinateMode = defaultCoordinateMode;
            this.coordinateRanges = defaultCoordinateRanges;

            this.selectTree = selectTree;
            this.selectedLocationIds = initialLocationIds; // TODO: select NYC by default
            this.caseData = aggCaseDataList;
            this.changingPositions = changingPositions;
            this.caseDataAggGroup = caseDataAggGroup;
            this.dateRange = initialDateRange;
            this.selectedRows = selectedRows;

            console.log('DATA INIT FINISHED');

            uiStoreInstance.onAggCaseDataFinished();
            uiStoreInstance.onCaseDataStateFinished();
          }
        );
      }
    );
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
    console.log('CHANGE COORDINATE MODE', coordinateMode);
    console.log('SELECTED GENE:', selectedGene);
    console.log('SELECTED PROTEIN:', selectedProtein);
    console.log('SELECTED PRIMERS:', selectedPrimers);
    console.log('CUSTOM COORDINATES:', customCoordinates);

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
      },
      ({ aggCaseDataList, selectedRows }) => {
        this.caseData = aggCaseDataList;
        this.selectedRows = selectedRows;
        console.log('CASE_DATA FINISHED');

        this.updateAggCaseDataByGroup((suppressUIUpdate = false));
      }
    );
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

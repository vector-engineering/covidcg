import { observable, action, toJS } from 'mobx';

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
  changeCoordinateMode({ coordinateMode, selectedGene, selectedProtein }) {
    console.log('CHANGE COORDINATE MODE', coordinateMode);
    console.log('SELECTED GENE:', selectedGene);
    console.log('SELECTED PROTEIN:', selectedProtein);

    this.coordinateMode = coordinateMode;
    this.selectedGene = getGene(selectedGene);
    this.selectedProtein = getProtein(selectedProtein);

    // Set the coordinate range based off the coordinate mode
    if (coordinateMode === 'gene') {
      this.coordinateRanges = [
        [this.selectedGene.start, this.selectedGene.end],
      ];
    } else if (coordinateMode === 'protein') {
      this.coordinateRanges = this.selectedProtein.ranges;
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

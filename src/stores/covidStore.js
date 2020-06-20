import { observable, action, toJS } from 'mobx';

import {
  processCaseData,
  aggCaseDataByGroup,
} from '../utils/caseDataWorkerWrapper';
import {
  downloadAcknowledgements,
  downloadAggCaseData,
  downloadSequencesAndMetadata,
} from '../utils/downloadWorkerWrapper';
import { getGene, loadGeneOptions } from '../utils/gene';
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
  @observable genes = [];
  @observable selectedGene = {};
  @observable startPos = null;
  @observable endPos = null;
  @observable selectTree = [];
  @observable selectedLocationIds = [];
  @observable caseData = [];
  @observable changingPositions = {};
  @observable caseDataAggGroup = [];
  @observable dateRange = [];
  @observable selectedRows = [];

  constructor() {
    // Select the Spike gene by default
    let defaultGene = getGene('S');

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
        selectedGene: toJS(defaultGene),
        groupKey: toJS(initialGroupKey),
        dnaOrAa: toJS(initialDnaOrAa),
      },
      ({ aggCaseDataList, selectedRows }) => {
        uiStoreInstance.onAggCaseDataStarted();

        aggCaseDataByGroup(
          {
            caseData: toJS(aggCaseDataList),
            selectedGene: toJS(defaultGene),
            groupKey: toJS(initialGroupKey),
            dnaOrAa: toJS(initialDnaOrAa),
            dateRange: toJS(initialDateRange),
          },
          ({ changingPositions, caseDataAggGroup }) => {
            this.groupKey = initialGroupKey;
            this.dnaOrAa = initialDnaOrAa;
            this.genes = loadGeneOptions();
            this.selectedGene = defaultGene;
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
  selectGene(_selectedGene) {
    // im not sure what this var actually is
    console.log('SELECT_GENE', _selectedGene);
    this.selectedGene = getGene(_selectedGene);

    // Get matching clade_ids
    //let lineages = getLineagesFromGene(this.selectedGene);
    //this.selectedLineages = lineages;

    this.updateCaseData();
  }

  @action
  selectLocations(selectedNodes) {
    console.log('select locations: ', selectedNodes);
    this.selectedLocationIds = getLocationIds(selectedNodes);

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
        selectedGene: toJS(this.selectedGene),
        groupKey: toJS(this.groupKey),
        dnaOrAa: toJS(this.dnaOrAa),
        dateRange: toJS(this.dateRange),
      },
      ({ caseDataAggGroup, changingPositions }) => {
        // console.log(caseDataAggGroup);
        this.caseDataAggGroup = caseDataAggGroup;
        this.changingPositions = changingPositions;
        console.log('AGG_CASE_DATA FINISHED');

        suppressUIUpdate ? null : uiStoreInstance.onAggCaseDataFinished();
        suppressUIUpdate ? null : uiStoreInstance.onCaseDataStateFinished();
      }
    );
  }

  @action
  updateCaseData(suppressUIUpdate = false) {
    suppressUIUpdate ? null : uiStoreInstance.onCaseDataStateStarted();

    processCaseData(
      {
        selectedLocationIds: toJS(this.selectedLocationIds),
        selectedGene: toJS(this.selectedGene),
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

  @action
  downloadSequencesAndMetadata() {
    downloadSequencesAndMetadata(
      { selectedRows: toJS(this.selectedRows) },
      (res) => {
        //console.log(res);
        downloadBlobURL(
          res.blobURL,
          generateSelectionString(
            'sequences',
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

import { observable, action, toJS } from 'mobx';

import {
  processCaseData,
  aggCaseDataByGroup,
} from '../utils/caseDataWorkerWrapper';
import { getGene, loadGeneOptions } from '../utils/gene';
//import { getLineagesFromGene } from '../utils/lineageData';
import {
  loadSelectTree,
  getLocationByNameAndLevel,
  getLocationIds,
} from '../utils/location';
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

    processCaseData(
      {
        selectedLocationIds: toJS(initialLocationIds),
        selectedGene: toJS(defaultGene),
        groupKey: toJS(initialGroupKey),
        dnaOrAa: toJS(initialDnaOrAa),
      },
      (caseData) => {
        uiStoreInstance.onDataChangeFinished();
        this.updateAggCaseDataByGroup();

        aggCaseDataByGroup(
          {
            caseData: toJS(caseData),
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
            this.caseData = caseData;
            this.changingPositions = changingPositions;
            this.caseDataAggGroup = caseDataAggGroup;
            this.dateRange = initialDateRange;

            uiStoreInstance.onDataChangeFinished();
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
    console.log(selectedNodes);

    this.selectedLocationIds = getLocationIds(selectedNodes);

    this.updateCaseData();
  }

  @action
  selectDateRange(_dateRange) {
    this.dateRange = _dateRange;
    this.updateAggCaseDataByGroup();
  }

  @action
  updateAggCaseDataByGroup() {
    uiStoreInstance.onDataChangeStarted();
    aggCaseDataByGroup(
      {
        caseData: toJS(this.caseData),
        selectedGene: toJS(this.selectedGene),
        groupKey: toJS(this.groupKey),
        dnaOrAa: toJS(this.dnaOrAa),
        dateRange: toJS(this.dateRange),
      },
      ({ caseDataAggGroup, changingPositions }) => {
        this.caseDataAggGroup = caseDataAggGroup;
        this.changingPositions = changingPositions;
        uiStoreInstance.onDataChangeFinished();
      }
    );
  }

  @action
  updateCaseData() {
    uiStoreInstance.onDataChangeStarted();

    processCaseData(
      {
        selectedLocationIds: toJS(this.selectedLocationIds),
        selectedGene: toJS(this.selectedGene),
        groupKey: toJS(this.groupKey),
        dnaOrAa: toJS(this.dnaOrAa),
      },
      (data) => {
        this.caseData = data;
        console.log('SELECT_LOCATIONS FINISHED');

        uiStoreInstance.onDataChangeFinished();
        this.updateAggCaseDataByGroup();
      }
    );
  }
}

export default ObservableCovidStore;

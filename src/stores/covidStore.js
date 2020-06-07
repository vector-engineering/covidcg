import { observable, action } from 'mobx';

import { processCaseData, aggCaseDataByGroup } from '../utils/caseData';
import { getGene, loadGeneOptions } from '../utils/gene';
//import { getLineagesFromGene } from '../utils/lineageData';
import {
  loadSelectTree,
  getLocationByNameAndLevel,
  getLocationIds,
} from '../utils/location';

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

    // Load case data
    let caseData = processCaseData(
      initialLocationIds,
      defaultGene,
      initialGroupKey,
      initialDnaOrAa
    );
    let { caseDataAggGroup, changingPositions } = aggCaseDataByGroup(
      caseData,
      defaultGene,
      initialGroupKey,
      initialDnaOrAa,
      initialDateRange
    );

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
    console.log('SELECT_LOCATIONS');

    this.selectedLocationIds = getLocationIds(selectedNodes);

    console.log('Location IDs:', this.selectedLocationIds);

    this.updateCaseData();
  }

  @action
  selectDateRange(_dateRange) {
    this.dateRange = _dateRange;
    this.updateAggCaseDataByGroup();
  }

  @action
  updateAggCaseDataByGroup() {
    let { caseDataAggGroup, changingPositions } = aggCaseDataByGroup(
      this.caseData,
      this.selectedGene,
      this.groupKey,
      this.dnaOrAa,
      this.dateRange
    );

    this.caseDataAggGroup = caseDataAggGroup;
    this.changingPositions = changingPositions;
  }

  @action
  updateCaseData() {
    this.caseData = processCaseData(
      this.selectedLocationIds,
      this.selectedGene,
      this.groupKey,
      this.dnaOrAa
    );
    this.updateAggCaseDataByGroup();
  }
}

export default ObservableCovidStore;

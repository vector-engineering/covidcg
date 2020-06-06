import { observable, computed, action } from 'mobx';
import React from 'react';
import { useLocalStore } from 'mobx-react'; // 6.x or mobx-react-lite@1.4.0

import _ from 'underscore';

import { getGene, loadGeneOptions } from '../utils/gene';

import {
  loadSelectTree,
  getLocationByNameAndLevel,
  getLocationIds,
} from '../utils/location';

import {
  loadCaseData,
  processCaseData,
  aggCaseDataByLineage,
} from '../utils/caseData';

import { loadLineageData, getLineagesFromGene } from '../utils/lineageData';

class ObservableCovidStore {
  @observable genes = [];
  @observable selectedGene = null;
  @observable startPos = null;
  @observable endPos = null;
  @observable selectTree = [];
  @observable selectedLocationIds = []; // TODO: select NYC by default
  @observable selectedLineages = [];
  @observable initialLineageData = {};
  @observable initialCaseData = {};
  @observable caseData = {};
  @observable changingPositions = [];
  @observable caseDataAggLineageList = [];
  @observable dateRange = [];

  constructor() {
    // Load data
    let initialCaseData = loadCaseData();
    let initialLineageData = loadLineageData();

    // Process case data - turn date strings into date objects
    let processedCaseData = _.map(initialCaseData, (row) => {
      row.date = new Date(row.date).getTime();
      return row;
    });

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

    let initialLocationIds = NYCLocationId;
    let initialLineages = getLineagesFromGene(defaultGene);

    // No initial date range
    let initialDateRange = [-1, -1];

    // Load case data
    let caseData = processCaseData(processedCaseData, initialLocationIds);
    let { caseDataAggLineageList, changingPositions } = aggCaseDataByLineage(
      caseData,
      initialLineageData,
      defaultGene.start,
      defaultGene.end,
      initialDateRange
    );

    this.genes = loadGeneOptions();
    this.selectedGene = defaultGene.gene;
    this.startPos = defaultGene.start;
    this.endPos = defaultGene.end;
    this.selectTree = selectTree;
    this.selectedLocationIds = initialLocationIds; // TODO: select NYC by default
    this.selectedLineages = initialLineages;
    this.initialLineageData = initialLineageData;
    this.initialCaseData = processedCaseData;
    this.caseData = caseData;
    this.changingPositions = changingPositions;
    this.caseDataAggLineageList = caseDataAggLineageList;
    this.dateRange = initialDateRange;
  }

  @action
  selectGene(_selectedGene) {
    // im not sure what this var actually is
    console.log('SELECT_GENE', _selectedGene);

    let selectedGene = getGene(_selectedGene);

    this.startPos = selectedGene.start;
    this.endPos = selectedGene.end;

    // Get matching clade_ids
    let lineages = getLineagesFromGene(selectedGene);
    this.selectedLineages = lineages;

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
    let { caseDataAggLineageList, changingPositions } = aggCaseDataByLineage(
      this.caseData,
      this.initialLineageData,
      this.startPos,
      this.endPos,
      this.dateRange
    );

    this.caseDataAggLineageList = caseDataAggLineageList;
    this.changingPositions = changingPositions;
  }

  @action
  updateCaseData() {
    this.caseData = processCaseData(
      this.initialCaseData,
      this.selectedLocationIds
    );
  }
}

export const CovidStoreContext = React.createContext(null);

// eslint-disable-next-line react/prop-types
export const CovidStoreProvider = ({ children }) => {
  const store = useLocalStore(() => {
    return new ObservableCovidStore();
  });
  return (
    <CovidStoreContext.Provider value={store}>
      {children}
    </CovidStoreContext.Provider>
  );
};

export const useCovidStore = () => {
  const store = React.useContext(CovidStoreContext);
  if (!store) {
    throw new Error('useStore must be used within a StoreProvider.');
  }
  return store;
};

export const connectCovidStore = (Component) => {
  // eslint-disable-next-line react/display-name
  return (props) => {
    const covidStore = useCovidStore();
    return <Component covidStore={covidStore} {...props} />;
  };
};

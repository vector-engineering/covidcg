import _ from 'underscore';

import {
  getGene,
  loadGeneOptions
} from '../utils/gene';

import {
  loadSelectTree,
  getLocationByNameAndLevel,
  getLocationIds
} from '../utils/location';

import {
  loadCaseData,
  processCaseData,
  aggCaseDataByLineage
} from '../utils/caseData';

import {
  loadLineageData,
  getLineagesFromGene
} from '../utils/lineageData'

// Load data
let initialCaseData = loadCaseData();
let initialLineageData = loadLineageData();

// Process case data - turn date strings into date objects
let processedCaseData = _.map(initialCaseData, row => {
  row.date = new Date(row.date).getTime();
  return row;
});

// Select the Spike gene by default
let defaultGene = getGene('S');

let selectTree = loadSelectTree()
// Select NYC by default
let NYCNode = getLocationByNameAndLevel(selectTree, 'New York City', 'location');
NYCNode[0].checked = true;
let NYCLocationId = getLocationIds(NYCNode);

let initialLocationIds = NYCLocationId;
let initialLineages = getLineagesFromGene(defaultGene);

// No initial date range
let initialDateRange = [-1, -1];

// Load case data
let caseData = processCaseData(processedCaseData, initialLocationIds);
let { caseDataAggLineageList, changingPositions } = aggCaseDataByLineage(caseData, initialLineageData, defaultGene.start, defaultGene.end, initialDateRange);



export default {
  covid: {
    genes: loadGeneOptions(),
    selectedGene: defaultGene.gene,
    startPos: defaultGene.start,
    endPos: defaultGene.end,
    selectTree: selectTree,
    selectedLocationIds: initialLocationIds, // TODO: select NYC by default
    selectedLineages: initialLineages,
    initialLineageData: initialLineageData,
    initialCaseData: processedCaseData,
    caseData: caseData,
    changingPositions: changingPositions,
    caseDataAggLineageList: caseDataAggLineageList,
    dateRange: initialDateRange
  }
};

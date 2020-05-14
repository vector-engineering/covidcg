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
  aggCaseDataByClade
} from '../utils/caseData';

import {
  getCladesFromGene, loadCladeData
} from '../utils/cladeData'

// Load data
let initialCaseData = loadCaseData();
let initialCladeData = loadCladeData();

// Process case data - turn date strings into date objects
let processedCaseData = _.map(initialCaseData, row => {
  row.date = Date.parse(row.date);
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
let initialClades = getCladesFromGene(defaultGene);
//let initialCladeIds = _.map(initialClades, clade => clade.index);

// No initial date range
let initialDateRange = [-1, -1];

// Load case data
let caseData = processCaseData(processedCaseData, initialCladeData, initialLocationIds);
let { caseDataAggCladeList, changingPositions } = aggCaseDataByClade(caseData, initialCladeData, defaultGene.start, defaultGene.end, initialDateRange);



export default {
  covid: {
    genes: loadGeneOptions(),
    selectedGene: defaultGene.gene,
    startPos: defaultGene.start,
    endPos: defaultGene.end,
    selectTree: selectTree,
    selectedLocationIds: initialLocationIds, // TODO: select NYC by default
    selectedClades: initialClades,
    initialCladeData: initialCladeData,
    initialCaseData: processedCaseData,
    caseData: caseData,
    changingPositions: changingPositions,
    caseDataAggCladeList: caseDataAggCladeList,
    dateRange: initialDateRange
  }
};

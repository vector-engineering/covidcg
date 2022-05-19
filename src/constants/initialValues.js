import { getGene, getProtein } from '../utils/gene_protein';
import { intToISO, ISOToInt } from '../utils/date';
import { config } from '../config';

import {
  GROUP_MUTATION,
  DNA_OR_AA,
  COORDINATE_MODES,
  NORM_MODES,
  COUNT_MODES,
  DATE_BINS,
  SORT_DIRECTIONS,
  TREE_COLOR_MODES,
  LOW_FREQ_FILTER_TYPES,
} from './defs.json';

const today = intToISO(new Date().getTime());
const lastNDays = 30;

let _configStore, _groupDataStore, _plotSettingsStore;

if (config['virus'] === 'sars2') {
  _configStore = {
    groupKey: GROUP_MUTATION,
    dnaOrAa: DNA_OR_AA.AA,

    // Select the Spike gene and nsp13 protein by default
    selectedGene: getGene('S'),
    selectedProtein: getProtein('nsp12 - RdRp'),
    selectedPrimers: [],
    customCoordinates: [[8000, 12000]],
    customSequences: ['GACCCCAAAATCAGCGAAAT'],
    residueCoordinates: [[1, getGene('S').len_aa]],

    // Selecting the gene as the coordinate range by default
    coordinateMode: COORDINATE_MODES.COORD_GENE,

    // days * (24 hours/day) * (60 min/hour) * (60 s/min) * (1000 ms/s)
    startDate: intToISO(ISOToInt(today) - lastNDays * 24 * 60 * 60 * 1000),
    endDate: today,

    submStartDate: '',
    submEndDate: '',

    selectedLocationNodes: [],

    hoverGroup: null,
    selectedGroups: [],

    // Metadata filtering
    selectedMetadataFields: {},
    ageRange: [null, null],

    // Location tab
    hoverLocation: null,
    focusedLocations: [],
  };

  _plotSettingsStore = {
    groupStackLowFreqFilter: LOW_FREQ_FILTER_TYPES.GROUP_COUNTS,
    groupStackLowFreqValue: 20,
    groupStackNormMode: NORM_MODES.NORM_COUNTS,
    groupStackCountMode: COUNT_MODES.COUNT_NEW,
    groupStackDateBin: DATE_BINS.DATE_BIN_DAY,

    locationDateNormMode: NORM_MODES.NORM_PERCENTAGES,
    locationDateCountMode: COUNT_MODES.COUNT_CUMULATIVE,
    locationDateDateBin: DATE_BINS.DATE_BIN_DAY,

    locationGroupHideReference: true,

    cooccurrenceNormMode: NORM_MODES.NORM_COUNTS,

    // SURVEILLANCE PLOT
    surveillanceMode: 'lineage',
    surveillanceSortField: 'counts', // 'group' or 'counts'
    surveillanceSortDirection: SORT_DIRECTIONS.SORT_DESC,
    surveillanceDisplayMinCounts: 5,
    surveillanceDisplayMinPercent: 0.01,
    surveillanceSigMinCounts: 10,
    surveillanceSigMinPercent: 0.02,
    surveillanceSigMinR: 0.3,
    surveillanceLegendHover: [],
    surveillanceShowWarning: false,
    surveillanceShowSettings: false,

    // GROUP REPORT TAB
    reportTreeColorMode: TREE_COLOR_MODES.COLOR_LATEST,
    reportConsensusThreshold: 0.7,
    reportMutationListHideEmpty: true,
    reportMutationListHidden: ['ORF1a'], // By default, hide ORF1a
    reportStructureActiveProtein: 'S',
    reportStructurePdbId: '6ZGG',
    reportStructureActiveGroup: 'B.1.617.2',
  };

  _groupDataStore = {
    activeGroupType: Object.keys(config['group_cols'])[0],
    selectedGroups: ['BA.1', 'AY.4', 'B.1.617.2', 'B.1.1.7', 'B.1.351', 'P.2'],
    groupMutationType: 'protein_aa',
    activeReference: '',
  };
} else if (config['virus'] === 'rsv') {
  const startingReference = 'NC_038235.1';

  _configStore = {
    groupKey: GROUP_MUTATION,
    dnaOrAa: DNA_OR_AA.AA,

    // Select the F gene and protein by default
    selectedGene: getGene('F', startingReference),
    selectedProtein: getProtein('F', startingReference),
    selectedPrimers: [],
    selectedReference: startingReference,
    customCoordinates: [[5648, 7550]],
    customSequences: ['GGTGTTGGATCTGCAATCGC'],
    residueCoordinates: [[1, getGene('F', startingReference).len_aa]],

    // Selecting the gene as the coordinate range by default
    coordinateMode: COORDINATE_MODES.COORD_GENE,

    // RSV displays all data by default
    startDate: config.min_date,
    endDate: today,

    submStartDate: '',
    submEndDate: '',

    selectedGroupFields: {
      subtype: ['A'],
    },
    selectedLocationNodes: [],

    hoverGroup: null,
    selectedGroups: [],

    // Metadata filtering
    selectedMetadataFields: {},
    ageRange: [null, null],

    // Location tab
    hoverLocation: null,
    focusedLocations: [],
  };
  _plotSettingsStore = {
    groupStackLowFreqFilter: LOW_FREQ_FILTER_TYPES.GROUP_COUNTS,
    groupStackLowFreqValue: 20,
    groupStackNormMode: NORM_MODES.NORM_COUNTS,
    groupStackCountMode: COUNT_MODES.COUNT_NEW,
    groupStackDateBin: DATE_BINS.DATE_BIN_MONTH,

    locationDateNormMode: NORM_MODES.NORM_PERCENTAGES,
    locationDateCountMode: COUNT_MODES.COUNT_CUMULATIVE,
    locationDateDateBin: DATE_BINS.DATE_BIN_MONTH,

    locationGroupHideReference: true,

    cooccurrenceNormMode: NORM_MODES.NORM_COUNTS,

    // SURVEILLANCE PLOT
    surveillanceMode: 'genotype',
    surveillanceSortField: 'counts', // 'group' or 'counts'
    surveillanceSortDirection: SORT_DIRECTIONS.SORT_DESC,
    surveillanceDisplayMinCounts: 5,
    surveillanceDisplayMinPercent: 0.01,
    surveillanceSigMinCounts: 10,
    surveillanceSigMinPercent: 0.02,
    surveillanceSigMinR: 0.15,
    surveillanceLegendHover: [],
    surveillanceShowWarning: false,
    surveillanceShowSettings: false,

    // GROUP REPORT TAB
    reportTreeColorMode: TREE_COLOR_MODES.COLOR_LATEST,
    reportConsensusThreshold: 0.7,
    reportMutationListHideEmpty: true,
    reportMutationListHidden: [], // By default, hide none
    reportStructureActiveProtein: 'F',
    reportStructurePdbId: '5UDE',
    reportStructureActiveGroup: 'ON1',
  };
  _groupDataStore = {
    activeGroupType: Object.keys(config['group_cols'])[0],
    selectedGroups: ['ON1'],
    groupMutationType: 'protein_aa',
    activeReference: 'NC_038235.1',
  };
}

export const configStore = _configStore;
export const plotSettingsStore = _plotSettingsStore;
export const groupDataStore = _groupDataStore;

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
  LITEMOL_STYLES,
} from './defs.json';

const today = intToISO(new Date().getTime());
const lastNDays = 30;

let _configStore, _groupDataStore, _plotSettingsStore;

if (config['virus'] === 'sars2') {
  const startingReference = '...';

  _configStore = {
    groupKey: GROUP_MUTATION,
    dnaOrAa: DNA_OR_AA.AA,

    // Select the Spike gene and nsp13 protein by default
    selectedGene: getGene(config['default_gene'], startingReference),
    selectedProtein: getProtein(config['default_protein'], startingReference),
    selectedPrimers: [],
    customCoordinates: [[8000, 12000]],
    customSequences: ['GACCCCAAAATCAGCGAAAT'],
    residueCoordinates: [
      [1, getGene(config['default_gene'], startingReference).len_aa],
    ],

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
    // LEGEND
    legendAdjustPartialSequences: false,

    // ENTROPY PLOT
    entropyYMode: NORM_MODES.NORM_COUNTS,
    entropyYPow: 0.5,

    // GROUP STACK PLOT
    groupStackLowFreqFilter: LOW_FREQ_FILTER_TYPES.GROUP_COUNTS,
    groupStackLowFreqValue: 20,
    groupStackNormMode: NORM_MODES.NORM_COUNTS,
    groupStackCountMode: COUNT_MODES.COUNT_NEW,
    groupStackDateBin: DATE_BINS.DATE_BIN_DAY,

    // MUTATION STRUCTURE VIEWER
    mutationStructurePdbId: '6ZGG',
    mutationStructureProteinStyle: LITEMOL_STYLES.SURFACE,
    mutationStructureNormMode: NORM_MODES.NORM_COVERAGE_ADJUSTED,

    // LOCATION DATE PLOT
    locationDateNormMode: NORM_MODES.NORM_PERCENTAGES,
    locationDateCountMode: COUNT_MODES.COUNT_CUMULATIVE,
    locationDateDateBin: DATE_BINS.DATE_BIN_DAY,

    // LOCATION GROUP PLOT
    locationGroupHideReference: true,

    // COOCCURRENCE PLOT
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
    reportStructureActiveProtein: config['default_protein'],
    reportStructurePdbId: '6ZGG',
    reportStructureActiveGroup: 'B.1.617.2',
    reportStructureProteinStyle: LITEMOL_STYLES.SURFACE,
  };

  _groupDataStore = {
    activeGroupType: Object.keys(config['group_cols'])[0],
    selectedGroups: ['BA.1', 'AY.4', 'B.1.617.2', 'B.1.1.7', 'B.1.351', 'P.2'],
    groupMutationType: 'protein_aa',
    activeReference: startingReference,
  };
} else if (config['virus'] === 'rsv') {
  const startingReference = 'NC_038235.1';

  _configStore = {
    groupKey: GROUP_MUTATION,
    dnaOrAa: DNA_OR_AA.AA,

    // Select the F gene and protein by default
    selectedGene: getGene(config['default_gene'], startingReference),
    selectedProtein: getProtein(config['default_protein'], startingReference),
    selectedPrimers: [],
    selectedReference: startingReference,
    customCoordinates: getGene(config['default_gene'], startingReference)
      .segments,
    customSequences: ['GGTGTTGGATCTGCAATCGC'],
    residueCoordinates: [
      [1, getGene(config['default_gene'], startingReference).len_aa],
    ],

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
    // LEGEND
    legendAdjustPartialSequences: true,

    // ENTROPY PLOT
    entropyYMode: NORM_MODES.NORM_COUNTS,
    entropyYPow: 0.5,

    // GROUP STACK PLOT
    groupStackLowFreqFilter: LOW_FREQ_FILTER_TYPES.GROUP_COUNTS,
    groupStackLowFreqValue: 50,
    groupStackNormMode: NORM_MODES.NORM_COUNTS,
    groupStackCountMode: COUNT_MODES.COUNT_NEW,
    groupStackDateBin: DATE_BINS.DATE_BIN_YEAR,

    // MUTATION STRUCTURE VIEWER
    // 5UDE, 3RRR
    mutationStructurePdbId: '5UDE',
    mutationStructureProteinStyle: LITEMOL_STYLES.SURFACE,
    mutationStructureNormMode: NORM_MODES.NORM_COVERAGE_ADJUSTED,

    // LOCATION DATE PLOT
    locationDateNormMode: NORM_MODES.NORM_PERCENTAGES,
    locationDateCountMode: COUNT_MODES.COUNT_CUMULATIVE,
    locationDateDateBin: DATE_BINS.DATE_BIN_YEAR,

    // LOCATION GROUP PLOT
    locationGroupHideReference: true,

    // COOCCURRENCE PLOT
    cooccurrenceNormMode: NORM_MODES.NORM_COUNTS,

    // SURVEILLANCE PLOT
    surveillanceMode: 'subtype',
    surveillanceSortField: 'counts', // 'group' or 'counts'
    surveillanceSortDirection: SORT_DIRECTIONS.SORT_DESC,
    surveillanceDisplayMinCounts: 0,
    surveillanceDisplayMinPercent: 0.0,
    surveillanceSigMinCounts: 0,
    surveillanceSigMinPercent: 0.0,
    surveillanceSigMinR: -1.0,
    surveillanceLegendHover: [],
    surveillanceShowWarning: false,
    surveillanceShowSettings: false,

    // GROUP REPORT TAB
    reportTreeColorMode: TREE_COLOR_MODES.COLOR_LATEST,
    reportConsensusThreshold: 0.7,
    reportMutationListHideEmpty: true,
    reportMutationListHidden: [], // By default, hide none
    reportStructureActiveProtein: config['default_protein'],
    reportStructurePdbId: '5UDE',
    reportStructureActiveGroup: 'ON1',
    reportStructureProteinStyle: LITEMOL_STYLES.SURFACE,
  };
  _groupDataStore = {
    activeGroupType: Object.keys(config['group_cols'])[0],
    selectedGroups: ['ON1'],
    groupMutationType: 'protein_aa',
    activeReference: startingReference,
  };
}

export const configStore = _configStore;
export const plotSettingsStore = _plotSettingsStore;
export const groupDataStore = _groupDataStore;

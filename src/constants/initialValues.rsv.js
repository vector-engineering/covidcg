import { getGene, getProtein } from '../utils/gene_protein';
import { intToISO } from '../utils/date';
import { config } from '../config';

import {
  GROUP_MUTATION,
  DNA_OR_AA,
  COORDINATE_MODES,
  MIN_DATE,
  NORM_MODES,
  COUNT_MODES,
  DATE_BINS,
  SORT_DIRECTIONS,
  TREE_COLOR_MODES,
  LOW_FREQ_FILTER_TYPES,
} from './defs.json';

const today = intToISO(new Date().getTime());

export default function values() {
  return {
    configStore: {
      groupKey: GROUP_MUTATION,
      dnaOrAa: DNA_OR_AA.AA,

      // Select the F gene and protein by default
      selectedGene: getGene('F', 'A'),
      selectedProtein: getProtein('F', 'A'),
      selectedPrimers: [],
      selectedReference: 'A',
      customCoordinates: [[5648, 7550]],
      customSequences: ['GGTGTTGGATCTGCAATCGC'],
      residueCoordinates: [[1, getGene('F', 'A').len_aa]],

      // Selecting the gene as the coordinate range by default
      coordinateMode: COORDINATE_MODES.COORD_GENE,

      // RSV displays all data by default
      startDate: MIN_DATE['RSV'],
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
    },
    plotSettingsStore: {
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
    },
    groupDataStore: {
      activeGroupType: Object.keys(config['group_cols'])[0],
      selectedGroups: ['ON1'],
      groupMutationType: 'protein_aa',
    },
  };
}
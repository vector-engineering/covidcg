import { getGene, getProtein } from '../utils/gene_protein';
import { intToISO, ISOToInt } from '../utils/date';

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
} from '../constants/defs.json';

const today = intToISO(new Date().getTime());
const lastNDays = 30; // By default, show only the last 1 month

export default function values() {
  return {
    configStore: {
      groupKey: GROUP_MUTATION,
      dnaOrAa: DNA_OR_AA.AA,

      // Select the F gene and protein by default
      selectedGene: getGene('F'),
      selectedProtein: getProtein('F'),
      selectedPrimers: [],
      selectedReference: 'A',
      customCoordinates: [[8000, 12000]],
      customSequences: ['GACCCCAAAATCAGCGAAAT'],
      residueCoordinates: [[1, getGene('F').len_aa]],

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
    },
    plotSettingsStore: {
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
      surveillanceMode: 'genotype',
      surveillanceSortField: 'counts', // 'group' or 'counts'
      surveillanceSortDirection: SORT_DIRECTIONS.SORT_DESC,
      surveillanceDisplayMinCounts: 5,
      surveillanceDisplayMinPercent: 0.01,
      surveillanceSigMinCounts: 10,
      surveillanceSigMinPercent: 0.02,
      surveillanceSigMinR: 0.3,
      surveillanceLegendHover: [],

      // GROUP REPORT TAB
      reportTreeColorMode: TREE_COLOR_MODES.COLOR_LATEST,
      reportConsensusThreshold: 0.7,
      reportMutationListHideEmpty: true,
      reportMutationListHidden: [], // By default, hide none
      reportStructureActiveProtein: 'F',
      reportStructurePdbId: '',
      reportStructureActiveGroup: 'A',
    },
  };
}

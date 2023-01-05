import { observable, action } from 'mobx';
import {
  NORM_MODES,
  COUNT_MODES,
  DATE_BINS,
  SORT_DIRECTIONS,
  TREE_COLOR_MODES,
  LOW_FREQ_FILTER_TYPES,
  LITEMOL_STYLES,
  SURV_GRAPH_MODE,
} from '../constants/defs.json';
import { plotSettingsStore as initialPlotSettingsStore } from '../constants/initialValues';

import { rootStoreInstance } from './rootStore';

import { updateURLFromParams } from '../utils/updateQueryParam';

export class PlotSettingsStore {
  initialValues = {};

  // LEGEND
  @observable legendAdjustPartialSequences = true;

  // ENTROPY PLOT
  @observable entropyYMode = NORM_MODES.NORM_COUNTS;
  @observable entropyYPow = 0.5;

  // GROUP STACK PLOT
  @observable groupStackLowFreqFilter = LOW_FREQ_FILTER_TYPES.GROUP_COUNTS;
  @observable groupStackLowFreqValue = 20;
  @observable groupStackNormMode = NORM_MODES.NORM_COUNTS;
  @observable groupStackCountMode = COUNT_MODES.COUNT_NEW;
  @observable groupStackDateBin = DATE_BINS.DATE_BIN_DAY;

  // MUTATION STRUCTURE VIEWER
  @observable mutationStructurePdbId = '';
  @observable mutationStructureProteinStyle = LITEMOL_STYLES.SURFACE;
  @observable mutationStructureNormMode = NORM_MODES.NORM_COVERAGE_ADJUSTED;

  mutationStructureAssemblies = [];
  mutationStructureActiveAssembly = '';
  mutationStructureEntities = [];

  // LOCATION DATE PLOT
  @observable locationDateNormMode = NORM_MODES.NORM_PERCENTAGES;
  @observable locationDateCountMode = COUNT_MODES.COUNT_CUMULATIVE;
  @observable locationDateDateBin = DATE_BINS.DATE_BIN_DAY;

  // LOCATION GROUP PLOT
  @observable locationGroupHideReference = true;

  // COOCCURRENCE PLOT
  @observable cooccurrenceNormMode = NORM_MODES.NORM_COUNTS;

  // SURVEILLANCE PLOT
  @observable surveillanceMode = '';
  @observable surveillanceGraphMode = SURV_GRAPH_MODE.LINE;
  @observable surveillanceShowWarning = true;
  @observable surveillanceShowSettings = false;
  @observable surveillanceSortField = '';
  @observable surveillanceSortDirection = SORT_DIRECTIONS.SORT_DESC;
  @observable surveillanceDisplayMinCounts = 5;
  @observable surveillanceDisplayMinPercent = 0.01;
  @observable surveillanceSigMinCounts = 5;
  @observable surveillanceSigMinPercent = 0.01;
  @observable surveillanceSigMinR = 0.3;
  @observable surveillanceLegendHover = [];

  // GROUP REPORT TAB
  @observable reportTreeColorMode = TREE_COLOR_MODES.COLOR_LATEST;
  @observable reportConsensusThreshold = 0.7;
  @observable reportMutationListHideEmpty = true;
  @observable reportMutationListHidden = [];
  @observable reportStructureActiveProtein = '';
  @observable reportStructurePdbId = '';
  // Actively selected group for the structural viewer
  @observable reportStructureActiveGroup = '';
  @observable reportStructureProteinStyle = LITEMOL_STYLES.SURFACE;

  reportStructureAssemblies = [];
  reportStructureActiveAssembly = '';
  reportStructureEntities = [];

  ignoreURLProperties = [
    'mutationStructureAssemblies',
    'mutationStructureActiveAssembly',
    'mutationStructureEntities',
    'reportStructureAssemblies',
    'reportStructureActiveAssembly',
    'reportStructureEntities',
    'surveillanceLegendHover',
  ];

  init() {
    this.initialValues = initialPlotSettingsStore;

    Object.keys(this.initialValues).forEach((key) => {
      this[key] = this.initialValues[key];
    });
  }

  @action
  applyPendingChanges = (pending) => {
    // Overwrite any of our fields here with the pending ones
    Object.keys(pending).forEach((field) => {
      if (
        field === 'reportMutationListHidden' &&
        typeof pending[field] === String
      ) {
        this[field] = pending[field].split(',');
      } else {
        this[field] = pending[field];
      }
    });

    this.updateURL(pending);
  };

  /*
   * Serialize store state into URL params
   */
  updateURL = (pending) => {
    const urlParams = rootStoreInstance.urlMonitor.urlParams;

    Object.keys(pending).forEach((field) => {
      if (field === 'reportMutationListHidden') {
        urlParams.set(field, pending[field].join(','));
      } else {
        urlParams.set(field, String(pending[field]));
      }

      if (
        JSON.stringify(pending[field]) ===
        JSON.stringify(this.initialValues[field])
      ) {
        // Only display non-default fields in the url
        urlParams.delete(field);
      } else if (this.ignoreURLProperties.includes(field)) {
        urlParams.delete(field);
      }
    });

    // Update URL
    updateURLFromParams(urlParams);

    rootStoreInstance.urlMonitor.urlParams = urlParams;
  };

  @action
  resetValues(values) {
    Object.keys(this.initialValues).forEach((key) => {
      if (key in values) {
        this[key] = values[key];
      } else {
        this[key] = this.initialValues[key];
      }
    });
  }
}

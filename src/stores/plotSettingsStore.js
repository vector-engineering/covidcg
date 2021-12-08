import { observable, action } from 'mobx';
import {
  NORM_MODES,
  COUNT_MODES,
  DATE_BINS,
  SORT_DIRECTIONS,
  TREE_COLOR_MODES,
  LOW_FREQ_FILTER_TYPES,
  LITEMOL_STYLES,
} from '../constants/defs.json';

export const initialValues = {
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
  surveillanceShowWarning: true,
  surveillanceShowSettings: false,
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
  reportMutationListHidden: ['ORF1a'], // By default, hide ORF1a
  reportStructureActiveProtein: 'S',
  // reportStructureActiveProtein: 'nsp5 - 3CLp',
  reportStructurePdbId: '6ZGG',
  // reportStructurePdbId: '7RFW',
  reportStructureActiveGroup: 'B.1.1.529',
  // reportStructureActiveGroup: 'B.1.351',
  reportStructureProteinStyle: LITEMOL_STYLES.SURFACE,
};

export class PlotSettingsStore {
  init() {}

  @observable groupStackLowFreqFilter = initialValues.groupStackLowFreqFilter;
  @observable groupStackLowFreqValue = initialValues.groupStackLowFreqValue;
  @observable groupStackNormMode = initialValues.groupStackNormMode;
  @observable groupStackCountMode = initialValues.groupStackCountMode;
  @observable groupStackDateBin = initialValues.groupStackDateBin;

  @action
  setGroupStackLowFreqFilter = (filterType) => {
    this.groupStackLowFreqFilter = filterType;
  };
  @action
  setGroupStackLowFreqValue = (filterValue) => {
    this.groupStackLowFreqValue = filterValue;
  };
  @action
  setGroupStackNormMode = (mode) => {
    this.groupStackNormMode = mode;
  };
  @action
  setGroupStackCountMode = (mode) => {
    this.groupStackCountMode = mode;
  };
  @action
  setGroupStackDateBin = (dateBin) => {
    this.groupStackDateBin = dateBin;
  };

  @observable locationDateNormMode = initialValues.locationDateNormMode;
  @observable locationDateCountMode = initialValues.locationDateCountMode;
  @observable locationDateDateBin = initialValues.locationDateDateBin;

  @action
  setLocationDateNormMode = (mode) => {
    this.locationDateNormMode = mode;
  };
  @action
  setLocationDateCountMode = (mode) => {
    this.locationDateCountMode = mode;
  };
  @action
  setLocationDateDateBin = (dateBin) => {
    this.locationDateDateBin = dateBin;
  };

  @observable locationGroupHideReference =
    initialValues.locationGroupHideReference;

  @action
  setLocationGroupHideReference = (hide) => {
    this.locationGroupHideReference = hide;
  };

  @observable cooccurrenceNormMode = initialValues.cooccurrenceNormMode;

  @action
  setCooccurrenceNormMode = (mode) => {
    this.cooccurrenceNormMode = mode;
  };

  // -----------------
  // SURVEILLANCE PLOT
  // -----------------

  @observable surveillanceMode = initialValues.surveillanceMode;
  @observable surveillanceShowWarning = initialValues.surveillanceShowWarning;
  @observable surveillanceShowSettings = initialValues.surveillanceShowSettings;
  @observable surveillanceSortField = initialValues.surveillanceSortField;
  @observable surveillanceSortDirection =
    initialValues.surveillanceSortDirection;
  @observable surveillanceDisplayMinCounts =
    initialValues.surveillanceDisplayMinCounts;
  @observable surveillanceDisplayMinPercent =
    initialValues.surveillanceDisplayMinPercent;
  @observable surveillanceSigMinCounts = initialValues.surveillanceSigMinCounts;
  @observable surveillanceSigMinPercent =
    initialValues.surveillanceSigMinPercent;
  @observable surveillanceSigMinR = initialValues.surveillanceSigMinR;
  @observable surveillanceLegendHover = initialValues.surveillanceLegendHover;

  @action
  setSurveillanceMode = (mode) => {
    this.surveillanceMode = mode;
  };
  @action
  setSurveillanceShowWarning = (show) => {
    this.surveillanceShowWarning = show;
  };
  @action
  setSurveillanceShowSettings = (show) => {
    this.surveillanceShowSettings = show;
  };
  @action
  setSurveillanceSortField = (field) => {
    this.surveillanceSortField = field;
  };
  @action
  setSurveillanceSortDirection = (direction) => {
    this.surveillanceSortDirection = direction;
  };
  @action
  setSurveillanceDisplayMinCounts = (counts) => {
    this.surveillanceDisplayMinCounts = counts;
  };
  @action
  setSurveillanceDisplayMinPercent = (percent) => {
    this.surveillanceDisplayMinPercent = percent;
  };
  @action
  setSurveillanceSigMinCounts = (counts) => {
    this.surveillanceSigMinCounts = counts;
  };
  @action
  setSurveillanceSigMinPercent = (percent) => {
    this.surveillanceSigMinPercent = percent;
  };
  @action
  setSurveillanceSigMinR = (r) => {
    this.surveillanceSigMinR = r;
  };
  @action
  setSurveillanceLegendHover = (hover) => {
    this.surveillanceLegendHover = hover;
  };

  // ----------------
  // GROUP REPORT TAB
  // ----------------

  @observable reportTreeColorMode = initialValues.reportTreeColorMode;
  @observable reportConsensusThreshold = initialValues.reportConsensusThreshold;
  @observable reportMutationListHideEmpty =
    initialValues.reportMutationListHideEmpty;
  @observable reportMutationListHidden = initialValues.reportMutationListHidden;
  @observable reportStructureActiveProtein =
    initialValues.reportStructureActiveProtein;
  @observable reportStructurePdbId = initialValues.reportStructurePdbId;
  // Actively selected group for the structural viewer
  @observable reportStructureActiveGroup =
    initialValues.reportStructureActiveGroup;
  @observable reportStructureProteinStyle =
    initialValues.reportStructureProteinStyle;

  @action
  setReportTreeColorMode = (mode) => {
    this.reportTreeColorMode = mode;
  };
  @action
  setReportConsensusThreshold = (thresh) => {
    this.reportConsensusThreshold = thresh;
  };
  @action
  setReportMutationListHideEmpty = (hideEmpty) => {
    this.reportMutationListHideEmpty = hideEmpty;
  };
  @action
  setReportMutationListHidden = (hidden) => {
    this.reportMutationListHidden = hidden;
  };
  @action
  toggleReportMutationListHiddenItem = (itemName) => {
    let hidden = this.reportMutationListHidden;

    // If the item is not in the list yet, then add it
    if (hidden.indexOf(itemName) === -1) {
      hidden.push(itemName);
    } else {
      // If it does exist, then make a new array without this item
      hidden = hidden.filter((name) => name != itemName);
    }

    this.reportMutationListHidden = hidden;
  };
  @action
  setReportStructureActiveProtein = (proteinName) => {
    this.reportStructureActiveProtein = proteinName;
  };
  @action
  setReportStructurePdbId = (pdbId) => {
    this.reportStructurePdbId = pdbId;
  };
  @action
  setReportStructureActiveGroup = (group) => {
    this.reportStructureActiveGroup = group;
  };
  @action
  setReportStructureProteinStyle = (style) => {
    this.reportStructureProteinStyle = style;
  };

  @action
  resetValues(values) {
    Object.keys(initialValues).forEach((key) => {
      if (key in values) {
        this[key] = values[key];
      } else {
        this[key] = initialValues[key];
      }
    });
  }
}

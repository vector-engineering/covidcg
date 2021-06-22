import { observable, action } from 'mobx';
import {
  NORM_MODES,
  COUNT_MODES,
  DATE_BINS,
  COLOR_MODES,
  COMPARE_MODES,
  COMPARE_COLORS,
  SORT_DIRECTIONS,
} from '../constants/defs.json';

export const initialPlotSettingsValues = {
  groupStackNormMode: NORM_MODES.NORM_COUNTS,
  groupStackCountMode: COUNT_MODES.COUNT_NEW,
  groupStackDateBin: DATE_BINS.DATE_BIN_DAY,

  locationDateNormMode: NORM_MODES.NORM_PERCENTAGES,
  locationDateCountMode: COUNT_MODES.COUNT_CUMULATIVE,
  locationDateDateBin: DATE_BINS.DATE_BIN_DAY,

  locationGroupHideReference: true,

  tableColorMode: COLOR_MODES.COLOR_MODE_COMPARE,
  tableCompareMode: COMPARE_MODES.COMPARE_MODE_MISMATCH,
  tableCompareColor: COMPARE_COLORS.COLOR_MODE_ZAPPO,
  tableSortColumn: 'counts',
  tableSortDirection: SORT_DIRECTIONS.SORT_DESC,

  cooccurrenceNormMode: NORM_MODES.NORM_COUNTS,

  // SURVEILLANCE PLOT
  surveillanceMode: 'lineage',
  surveillanceSortField: 'group',
  surveillanceSortDirection: SORT_DIRECTIONS.SORT_ASC,
  surveillanceDisplayMinCounts: 5,
  surveillanceDisplayMinPercent: 0.01,
  surveillanceSigMinCounts: 10,
  surveillanceSigMinPercent: 0.05,
  surveillanceSigMinR: 0.5,
  surveillanceLegendHover: [],

  // GROUP REPORT TAB
  reportConsensusThreshold: 0.7,
  reportMutationListHidden: ['ORF1a'], // By default, hide Orf1a
};

export class PlotSettingsStore {
  @observable groupStackNormMode = initialPlotSettingsValues.groupStackNormMode;
  @observable groupStackCountMode =
    initialPlotSettingsValues.groupStackCountMode;
  @observable groupStackDateBin = initialPlotSettingsValues.groupStackDateBin;

  init() {}

  @action
  setGroupStackNormMode(mode) {
    this.groupStackNormMode = mode;
  }
  @action
  setGroupStackCountMode(mode) {
    this.groupStackCountMode = mode;
  }
  @action
  setGroupStackDateBin(dateBin) {
    this.groupStackDateBin = dateBin;
  }

  @observable locationDateNormMode =
    initialPlotSettingsValues.locationDateNormMode;
  @observable locationDateCountMode =
    initialPlotSettingsValues.locationDateCountMode;
  @observable locationDateDateBin =
    initialPlotSettingsValues.locationDateDateBin;

  @action
  setLocationDateNormMode(mode) {
    this.locationDateNormMode = mode;
  }
  @action
  setLocationDateCountMode(mode) {
    this.locationDateCountMode = mode;
  }
  @action
  setLocationDateDateBin(dateBin) {
    this.locationDateDateBin = dateBin;
  }

  @observable locationGroupHideReference =
    initialPlotSettingsValues.locationGroupHideReference;

  @action
  setLocationGroupHideReference(hide) {
    this.locationGroupHideReference = hide;
  }

  @observable tableColorMode = initialPlotSettingsValues.tableColorMode;
  @observable tableCompareMode = initialPlotSettingsValues.tableCompareMode;
  @observable tableCompareColor = initialPlotSettingsValues.tableCompareColor;
  @observable tableSortColumn = initialPlotSettingsValues.tableSortColumn;
  @observable tableSortDirection = initialPlotSettingsValues.tableSortDirection;

  @action
  setTableColorMode(mode) {
    this.tableColorMode = mode;
  }
  @action
  setTableCompareMode(mode) {
    this.tableCompareMode = mode;
  }
  @action
  setTableCompareColor(color) {
    this.tableCompareColor = color;
  }
  @action
  setTableSort(col, dir) {
    this.tableSortColumn = col;
    this.tableSortDirection = dir;
  }

  @observable cooccurrenceNormMode =
    initialPlotSettingsValues.cooccurrenceNormMode;

  @action
  setCooccurrenceNormMode(mode) {
    this.cooccurrenceNormMode = mode;
  }

  // -----------------
  // SURVEILLANCE PLOT
  // -----------------

  @observable surveillanceMode = initialPlotSettingsValues.surveillanceMode;
  @observable surveillanceSortField =
    initialPlotSettingsValues.surveillanceSortField;
  @observable surveillanceSortDirection =
    initialPlotSettingsValues.surveillanceSortDirection;
  @observable surveillanceDisplayMinCounts =
    initialPlotSettingsValues.surveillanceDisplayMinCounts;
  @observable surveillanceDisplayMinPercent =
    initialPlotSettingsValues.surveillanceDisplayMinPercent;
  @observable surveillanceSigMinCounts =
    initialPlotSettingsValues.surveillanceSigMinCounts;
  @observable surveillanceSigMinPercent =
    initialPlotSettingsValues.surveillanceSigMinPercent;
  @observable surveillanceSigMinR =
    initialPlotSettingsValues.surveillanceSigMinR;
  @observable surveillanceLegendHover =
    initialPlotSettingsValues.surveillanceLegendHover;

  @action
  setSurveillanceMode(mode) {
    this.surveillanceMode = mode;
  }
  @action
  setSurveillanceSortField(field) {
    this.surveillanceSortField = field;
  }
  @action
  setSurveillanceSortDirection(direction) {
    this.surveillanceSortDirection = direction;
  }
  @action
  setSurveillanceDisplayMinCounts(counts) {
    this.surveillanceDisplayMinCounts = counts;
  }
  @action
  setSurveillanceDisplayMinPercent(percent) {
    this.surveillanceDisplayMinPercent = percent;
  }
  @action
  setSurveillanceSigMinCounts(counts) {
    this.surveillanceSigMinCounts = counts;
  }
  @action
  setSurveillanceSigMinPercent(percent) {
    this.surveillanceSigMinPercent = percent;
  }
  @action
  setSurveillanceSigMinR(r) {
    this.surveillanceSigMinR = r;
  }
  @action
  setSurveillanceLegendHover(hover) {
    this.surveillanceLegendHover = hover;
  }

  // ----------------
  // GROUP REPORT TAB
  // ----------------

  @observable reportConsensusThreshold =
    initialPlotSettingsValues.reportConsensusThreshold;
  @observable reportMutationListHidden =
    initialPlotSettingsValues.reportMutationListHidden;

  @action
  setReportConsensusThreshold(thresh) {
    this.reportConsensusThreshold = thresh;
  }
  @action
  setReportMutationListHidden(hidden) {
    this.reportMutationListHidden = hidden;
  }
  @action
  toggleReportMutationListHiddenItem(itemName) {
    let hidden = this.reportMutationListHidden;

    // If the item is not in the list yet, then add it
    if (hidden.indexOf(itemName) === -1) {
      hidden.push(itemName);
    } else {
      // If it does exist, then make a new array without this item
      hidden = hidden.filter((name) => name != itemName);
    }

    this.reportMutationListHidden = hidden;
  }

  @action
  resetValues(values) {
    Object.keys(initialPlotSettingsValues).forEach((key) => {
      if (key in values) {
        this[key] = values[key];
      } else {
        this[key] = initialPlotSettingsValues[key];
      }
    });
  }
}

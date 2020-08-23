import { observable, action } from 'mobx';
import {
  NORM_MODES,
  COUNT_MODES,
  DATE_BINS,
  COLOR_MODES,
  COMPARE_MODES,
  COMPARE_COLORS,
  SORT_DIRECTIONS,
} from '../constants/plotSettings';

export const initialPlotSettingsValues = {
  groupStackNormMode: NORM_MODES.NORM_COUNTS,
  groupStackCountMode: COUNT_MODES.COUNT_NEW,
  groupStackDateBin: DATE_BINS.DATE_BIN_DAY,

  locationDateNormMode: NORM_MODES.NORM_PERCENTAGES,
  locationDateCountMode: COUNT_MODES.COUNT_CUMULATIVE,
  locationDateDateBin: DATE_BINS.DATE_BIN_DAY,

  tableColorMode: COLOR_MODES.COLOR_MODE_COMPARE,
  tableCompareMode: COMPARE_MODES.COMPARE_MODE_MISMATCH,
  tableCompareColor: COMPARE_COLORS.COLOR_MODE_ZAPPO,
  tableSortColumn: 'cases_sum',
  tableSortDirection: SORT_DIRECTIONS.SORT_DESC,

  cooccurrenceNormMode: NORM_MODES.NORM_COUNTS,
};

class ObservablePlotSettingsStore {
  @observable groupStackNormMode = initialPlotSettingsValues.groupStackNormMode;
  @observable groupStackCountMode =
    initialPlotSettingsValues.groupStackCountMode;
  @observable groupStackDateBin = initialPlotSettingsValues.groupStackDateBin;

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

export default ObservablePlotSettingsStore;

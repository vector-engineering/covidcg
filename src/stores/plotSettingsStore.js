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

class ObservablePlotSettingsStore {
  @observable groupStackNormMode = NORM_MODES.NORM_COUNTS;
  @observable groupStackCountMode = COUNT_MODES.COUNT_NEW;
  @observable groupStackDateBin = DATE_BINS.DATE_BIN_DAY;

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

  @observable locationDateNormMode = NORM_MODES.NORM_PERCENTAGES;
  @observable locationDateCountMode = COUNT_MODES.COUNT_CUMULATIVE;
  @observable locationDateDateBin = DATE_BINS.DATE_BIN_DAY;

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

  // Color by 'compare': Comparison to reference, or 'code': With a defined color code
  @observable tableColorMode = COLOR_MODES.COLOR_MODE_COMPARE;
  @observable tableCompareMode = COMPARE_MODES.COMPARE_MODE_MISMATCH;
  @observable tableCompareColor = COMPARE_COLORS.COMPARE_COLOR_YELLOW;
  @observable tableSortColumn = 'cases_sum';
  @observable tableSortDirection = SORT_DIRECTIONS.SORT_DESC;

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
}

export default ObservablePlotSettingsStore;

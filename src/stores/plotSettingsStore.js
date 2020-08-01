import { observable, action, toJS } from 'mobx';

const NORM_COUNTS = 'NORM_COUNTS';
const NORM_PERCENTAGES = 'NORM_PERCENTAGES';
export const NORM_MODES = {
  NORM_COUNTS,
  NORM_PERCENTAGES,
};

const COUNT_NEW = 'COUNT_NEW';
const COUNT_CUMULATIVE = 'COUNT_CUMULATIVE';
export const COUNT_MODES = {
  COUNT_NEW,
  COUNT_CUMULATIVE,
};

const DATE_BIN_DAY = 'DATE_BIN_DAY';
const DATE_BIN_WEEK = 'DATE_BIN_WEEK';
const DATE_BIN_MONTH = 'DATE_BIN_MONTH';
export const DATE_BINS = {
  DATE_BIN_DAY,
  DATE_BIN_WEEK,
  DATE_BIN_MONTH,
};

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
}

export default ObservablePlotSettingsStore;

import { observable, action } from 'mobx';
import {
  NORM_MODES,
  COUNT_MODES,
  DATE_BINS,
  SORT_DIRECTIONS,
  TREE_COLOR_MODES,
  LOW_FREQ_FILTER_TYPES,
  // LITEMOL_STYLES,
} from '../constants/defs.json';
import { plotSettingsStore as initialPlotSettingsStore } from '../constants/initialValues';

export class PlotSettingsStore {
  initialValues = {};

  @observable groupStackLowFreqFilter = LOW_FREQ_FILTER_TYPES.GROUP_COUNTS;
  @observable groupStackLowFreqValue = 20;
  @observable groupStackNormMode = NORM_MODES.NORM_COUNTS;
  @observable groupStackCountMode = COUNT_MODES.COUNT_NEW;
  @observable groupStackDateBin = DATE_BINS.DATE_BIN_DAY;

  @observable locationDateNormMode = NORM_MODES.NORM_PERCENTAGES;
  @observable locationDateCountMode = COUNT_MODES.COUNT_CUMULATIVE;
  @observable locationDateDateBin = DATE_BINS.DATE_BIN_DAY;

  @observable locationGroupHideReference = true;

  @observable cooccurrenceNormMode = NORM_MODES.NORM_COUNTS;

  // -----------------
  // SURVEILLANCE PLOT
  // -----------------

  @observable surveillanceMode = '';
  @observable surveillanceSortField = '';
  @observable surveillanceSortDirection = SORT_DIRECTIONS.SORT_DESC;
  @observable surveillanceDisplayMinCounts = 5;
  @observable surveillanceDisplayMinPercent = 0.01;
  @observable surveillanceSigMinCounts = 5;
  @observable surveillanceSigMinPercent = 0.01;
  @observable surveillanceSigMinR = 0.3;
  @observable surveillanceLegendHover = [];

  @observable reportTreeColorMode = TREE_COLOR_MODES.COLOR_LATEST;
  @observable reportConsensusThreshold = 0.7;
  @observable reportMutationListHideEmpty = true;
  @observable reportMutationListHidden = [];
  @observable reportStructureActiveProtein = '';
  @observable reportStructurePdbId = '';
  // Actively selected group for the structural viewer
  @observable reportStructureActiveGroup = '';

  init() {
    this.initialValues = initialPlotSettingsStore;

    Object.keys(this.initialValues).forEach((key) => {
      this[key] = this.initialValues[key];
    });
  }

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

  @action
  setLocationGroupHideReference = (hide) => {
    this.locationGroupHideReference = hide;
  };

  @action
  setCooccurrenceNormMode = (mode) => {
    this.cooccurrenceNormMode = mode;
  };

  // -----------------
  // SURVEILLANCE PLOT
  // -----------------

  @observable surveillanceMode = this.initialValues.surveillanceMode;
  @observable surveillanceShowWarning =
    this.initialValues.surveillanceShowWarning;
  @observable surveillanceShowSettings =
    this.initialValues.surveillanceShowSettings;
  @observable surveillanceSortField = this.initialValues.surveillanceSortField;
  @observable surveillanceSortDirection =
    this.initialValues.surveillanceSortDirection;
  @observable surveillanceDisplayMinCounts =
    this.initialValues.surveillanceDisplayMinCounts;
  @observable surveillanceDisplayMinPercent =
    this.initialValues.surveillanceDisplayMinPercent;
  @observable surveillanceSigMinCounts =
    this.initialValues.surveillanceSigMinCounts;
  @observable surveillanceSigMinPercent =
    this.initialValues.surveillanceSigMinPercent;
  @observable surveillanceSigMinR = this.initialValues.surveillanceSigMinR;
  @observable surveillanceLegendHover =
    this.initialValues.surveillanceLegendHover;

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

  @observable reportTreeColorMode = this.initialValues.reportTreeColorMode;
  @observable reportConsensusThreshold =
    this.initialValues.reportConsensusThreshold;
  @observable reportMutationListHideEmpty =
    this.initialValues.reportMutationListHideEmpty;
  @observable reportMutationListHidden =
    this.initialValues.reportMutationListHidden;
  @observable reportStructureActiveProtein =
    this.initialValues.reportStructureActiveProtein;
  @observable reportStructurePdbId = this.initialValues.reportStructurePdbId;
  // Actively selected group for the structural viewer
  @observable reportStructureActiveGroup =
    this.initialValues.reportStructureActiveGroup;
  @observable reportStructureProteinStyle =
    this.initialValues.reportStructureProteinStyle;

  reportStructureAssemblies = this.initialValues.reportStructureAssemblies;
  reportStructureActiveAssembly =
    this.initialValues.reportStructureActiveAssembly;
  reportStructureEntities = this.initialValues.reportStructureEntities;

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
  setReportStructureAssemblies = (assemblies) => {
    this.reportStructureAssemblies = assemblies;
  };
  @action
  setReportStructureActiveAssembly = (assembly) => {
    this.reportStructureActiveAssembly = assembly;
  };
  @action
  setReportStructureEntities = (entities) => {
    this.reportStructureEntities = entities;
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

import { observable, action, toJS } from 'mobx';
import { ASYNC_STATES, TABS } from '../constants/defs.json';
import { rootStoreInstance } from './rootStore';
import { updateURLFromParams } from '../utils/updateQueryParam';

export const initialValues = {
  caseDataState: ASYNC_STATES.STARTED,
  mutationDataState: ASYNC_STATES.STARTED,
  cooccurrenceDataState: ASYNC_STATES.STARTED,
  downloadState: ASYNC_STATES.UNINITIALIZED,

  globalSequencingDataState: ASYNC_STATES.UNINITIALIZED,
  metadataFieldsState: ASYNC_STATES.UNINITIALIZED,

  reportGroupMutationFrequencyState: ASYNC_STATES.UNINITIALIZED,

  activeTab: TABS.TAB_EXAMPLE,
  keysPressed: [],
};

export class UIStore {
  @observable caseDataState = initialValues.caseDataState;
  @observable mutationDataState = initialValues.mutationDataState;
  @observable cooccurrenceDataState = initialValues.cooccurrenceDataState;
  @observable downloadState = initialValues.downloadState;

  // Flag for whether or not we have the latest set of metadata mappings
  // i.e., metadata key (integer) => metadata value (string)
  @observable metadataFieldState = initialValues.metadataFieldState;

  @observable reportGroupMutationFrequencyState =
    initialValues.reportGroupMutationFrequencyState;

  @observable activeTab = initialValues.activeTab;
  @observable keysPressed = initialValues.keysPressed;

  init() {}

  @action
  resetValues = (values) => {
    Object.keys(initialValues).forEach((key) => {
      if (key in values) {
        this[key] = values[key];
      } else {
        this[key] = initialValues[key];
      }
    });
  };

  @action
  onCaseDataStateStarted = () => {
    this.caseDataState = ASYNC_STATES.STARTED;
  };
  @action
  onCaseDataStateFinished = () => {
    this.caseDataState = ASYNC_STATES.SUCCEEDED;
  };
  @action
  onCaseDataStateErr = () => {
    this.caseDataState = ASYNC_STATES.FAILED;
  };

  @action
  onMutationDataStarted = () => {
    this.mutationDataState = ASYNC_STATES.STARTED;
  };
  @action
  onMutationDataFinished = () => {
    this.mutationDataState = ASYNC_STATES.SUCCEEDED;
  };
  @action
  onMutationDataErr = () => {
    this.mutationDataState = ASYNC_STATES.FAILED;
  };

  @action
  onCooccurrenceDataStarted = () => {
    this.cooccurrenceDataState = ASYNC_STATES.STARTED;
  };
  @action
  onCooccurrenceDataFinished = () => {
    this.cooccurrenceDataState = ASYNC_STATES.SUCCEEDED;
  };
  @action
  onCooccurrenceDataErr = () => {
    this.cooccurrenceDataState = ASYNC_STATES.FAILED;
  };

  @action
  onDownloadStarted = () => {
    this.downloadState = ASYNC_STATES.STARTED;
  };
  @action
  onDownloadFinished = () => {
    this.downloadState = ASYNC_STATES.SUCCEEDED;
  };
  @action
  onDownloadErr = () => {
    this.downloadState = ASYNC_STATES.FAILED;
  };

  @action
  onGlobalSequencingDataStarted = () => {
    this.globalSequencingDataState = ASYNC_STATES.STARTED;
  };
  @action
  onGlobalSequencingDataFinished = () => {
    this.globalSequencingDataState = ASYNC_STATES.SUCCEEDED;
  };
  @action
  onGlobalSequencingDataErr = () => {
    this.globalSequencingDataState = ASYNC_STATES.FAILED;
  };

  @action
  onMetadataFieldStarted = () => {
    this.metadataFieldState = ASYNC_STATES.STARTED;
  };
  @action
  onMetadataFieldFinished = () => {
    this.metadataFieldState = ASYNC_STATES.SUCCEEDED;
  };
  @action
  onMetadataFieldErr = () => {
    this.metadataFieldState = ASYNC_STATES.FAILED;
  };

  @action
  onReportGroupMutationFrequencyStarted = () => {
    this.reportGroupMutationFrequencyState = ASYNC_STATES.STARTED;
  };
  @action
  onReportGroupMutationFrequencyFinished = () => {
    this.reportGroupMutationFrequencyState = ASYNC_STATES.SUCCEEDED;
  };
  @action
  onReportGroupMutationFrequencyErr = () => {
    this.reportGroupMutationFrequencyState = ASYNC_STATES.FAILED;
  };

  @action
  setActiveTab = (tab) => {
    if (Object.values(TABS).includes(tab)) {
      this.activeTab = tab;
    } else {
      this.activeTab = TABS.TAB_EXAMPLE;
    }

    if (
      (this.activeTab === TABS.TAB_COMPARE_GROUPS ||
        this.activeTab === TABS.TAB_COMPARE_LOCATIONS) &&
      this.caseDataState === ASYNC_STATES.STARTED
    ) {
      rootStoreInstance.dataStore.fetchData();
    } else if (
      this.activeTab === TABS.TAB_GLOBAL_SEQUENCES &&
      (this.globalSequencingDataState !== ASYNC_STATES.SUCCEEDED ||
        this.globalSequencingDataState === ASYNC_STATES.STARTED)
    ) {
      rootStoreInstance.globalSequencingDataStore.fetchGlobalSequencingData();
    } else if (
      this.activeTab === TABS.TAB_GROUP_REPORT &&
      (this.reportGroupMutationFrequencyState !== ASYNC_STATES.SUCCEEDED ||
        this.reportGroupMutationFrequencyState === ASYNC_STATES.STARTED)
    ) {
      rootStoreInstance.groupDataStore.fetchGroupMutationFrequencyData({
        group: rootStoreInstance.groupDataStore.activeReportGroupType,
        mutationType: rootStoreInstance.groupDataStore.reportGroupMutationType,
        consensusThreshold: 0,
        selectedReference:
          rootStoreInstance.groupDataStore.reportActiveReference,
      });
    }

    rootStoreInstance.urlMonitor.urlParams.set('tab', this.activeTab);
    updateURLFromParams(rootStoreInstance.urlMonitor.urlParams);
  };

  @action
  setKeyPressed = (keyCode, value) => {
    let keysPressed = toJS(this.keysPressed);
    if (!keysPressed.includes(keyCode)) {
      if (value) {
        keysPressed.push(keyCode);
      }
    } else {
      if (!value) {
        keysPressed = keysPressed.filter((k) => k !== keyCode);
      }
    }
    this.keysPressed = keysPressed;
  };

  isKeyPressed(keyCode) {
    return this.keysPressed.includes(keyCode);
  }
}

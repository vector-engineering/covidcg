import { observable, action, toJS } from 'mobx';
import { ASYNC_STATES, TABS } from '../constants/defs.json';
import { rootStoreInstance } from './rootStore';
import { updateURLFromParams } from '../utils/updateQueryParam';

function removeItemAll(arr, value) {
  var i = 0;
  while (i < arr.length) {
    if (arr[i] === value) {
      arr.splice(i, 1);
    } else {
      ++i;
    }
  }
  return arr;
}

export const initialValues = {
  sidebarOpen: false,
  sidebarSelectedGroupKeys: [],

  caseDataState: ASYNC_STATES.STARTED,
  snvDataState: ASYNC_STATES.STARTED,
  cooccurrenceDataState: ASYNC_STATES.STARTED,
  downloadState: ASYNC_STATES.UNINITIALIZED,

  globalSequencingDataState: ASYNC_STATES.UNINITIALIZED,
  metadataFieldsState: ASYNC_STATES.UNINITIALIZED,

  groupSnvFrequencyState: ASYNC_STATES.UNINITIALIZED,

  activeTab: TABS.TAB_EXAMPLE,
  keysPressed: [],
};

export class UIStore {
  @observable sidebarOpen = initialValues.sidebarOpen;
  @observable sidebarSelectedGroupKeys = initialValues.sidebarSelectedGroupKeys;

  @observable caseDataState = initialValues.caseDataState;
  @observable snvDataState = initialValues.snvDataState;
  @observable cooccurrenceDataState = initialValues.cooccurrenceDataState;
  @observable downloadState = initialValues.downloadState;

  // Flag for whether or not we have the latest set of metadata mappings
  // i.e., metadata key (integer) => metadata value (string)
  @observable metadataFieldState = initialValues.metadataFieldState;

  @observable groupSnvFrequencyState = initialValues.groupSnvFrequencyState;

  @observable activeTab = initialValues.activeTab;
  @observable keysPressed = initialValues.keysPressed;

  init() {}

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
  onSnvDataStarted = () => {
    this.snvDataState = ASYNC_STATES.STARTED;
  };
  @action
  onSnvDataFinished = () => {
    this.snvDataState = ASYNC_STATES.SUCCEEDED;
  };
  @action
  onSnvDataErr = () => {
    this.snvDataState = ASYNC_STATES.FAILED;
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
  onGroupSnvFrequencyStarted = () => {
    this.groupSnvFrequencyState = ASYNC_STATES.STARTED;
  };
  @action
  onGroupSnvFrequencyFinished = () => {
    this.groupSnvFrequencyState = ASYNC_STATES.SUCCEEDED;
  };
  @action
  onGroupSnvFrequencyErr = () => {
    this.groupSnvFrequencyState = ASYNC_STATES.FAILED;
  };

  @action
  setSidebarOpen = () => {
    this.sidebarOpen = true;
  };

  @action
  setSidebarClosed = () => {
    this.sidebarOpen = false;
  };

  @action
  onSelectGroupForSidebar(groupKey) {
    this.sidebarSelectedGroupKeys.push(groupKey);
  }

  @action
  onRemoveGroupFromSidebar(groupKey) {
    this.sidebarSelectedGroupKeys = removeItemAll(
      this.sidebarSelectedGroupKeys,
      groupKey
    );
  }

  @action
  setActiveTab(tab) {
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
      (this.groupSnvFrequencyState !== ASYNC_STATES.SUCCEEDED ||
        this.groupSnvFrequencyState === ASYNC_STATES.STARTED)
    ) {
      rootStoreInstance.groupDataStore.fetchGroupSnvFrequencyData({
        group: rootStoreInstance.groupDataStore.activeGroupType,
        snvType: rootStoreInstance.groupDataStore.groupSnvType,
        consensusThreshold: 0,
      });
    }

    rootStoreInstance.configStore.urlParams.set('tab', this.activeTab);
    updateURLFromParams(rootStoreInstance.configStore.urlParams);
  }

  @action
  setKeyPressed(keyCode, value) {
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
  }

  isKeyPressed(keyCode) {
    return this.keysPressed.includes(keyCode);
  }
}

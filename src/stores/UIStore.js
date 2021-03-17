import { observable, action, toJS } from 'mobx';
import { ASYNC_STATES, TABS } from '../constants/defs.json';

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

export const initialUIValues = {
  sidebarOpen: false,
  sidebarSelectedGroupKeys: [],

  caseDataState: ASYNC_STATES.STARTED,
  snvDataState: ASYNC_STATES.STARTED,
  cooccurrenceDataState: ASYNC_STATES.STARTED,
  downloadState: ASYNC_STATES.SUCCEEDED,

  activeTab: TABS.TAB_EXAMPLE,
  keysPressed: [],
};

export class UIStore {
  @observable sidebarOpen = initialUIValues.sidebarOpen;
  @observable sidebarSelectedGroupKeys =
    initialUIValues.sidebarSelectedGroupKeys;

  @observable caseDataState = initialUIValues.caseDataState;
  @observable snvDataState = initialUIValues.snvDataState;
  @observable cooccurrenceDataState = initialUIValues.cooccurrenceDataState;
  @observable downloadState = initialUIValues.downloadState;

  @observable activeTab = initialUIValues.activeTab;
  @observable keysPressed = initialUIValues.keysPressed;

  init() {}

  @action
  resetValues(values) {
    Object.keys(initialUIValues).forEach((key) => {
      if (key in values) {
        this[key] = values[key];
      } else {
        this[key] = initialUIValues[key];
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

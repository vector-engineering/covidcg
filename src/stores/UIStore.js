import { observable, action, toJS } from 'mobx';
import { ASYNC_STATES, TABS } from '../constants/UI';

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
  aggCaseDataState: ASYNC_STATES.STARTED,
  snvDataState: ASYNC_STATES.STARTED,
  cooccurrenceDataState: ASYNC_STATES.STARTED,

  activeTab: TABS.TAB_EXAMPLE,
  keysPressed: [],
};

class ObservableUIStore {
  @observable sidebarOpen = initialUIValues.sidebarOpen;
  @observable sidebarSelectedGroupKeys =
    initialUIValues.sidebarSelectedGroupKeys;

  @observable caseDataState = initialUIValues.caseDataState;
  @observable aggCaseDataState = initialUIValues.aggCaseDataState;
  @observable snvDataState = initialUIValues.snvDataState;
  @observable cooccurrenceDataState = initialUIValues.cooccurrenceDataState;

  @observable activeTab = initialUIValues.activeTab;
  @observable keysPressed = initialUIValues.keysPressed;

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
  onAggCaseDataStarted = () => {
    this.aggCaseDataState = ASYNC_STATES.STARTED;
  };
  @action
  onAggCaseDataFinished = () => {
    this.aggCaseDataState = ASYNC_STATES.SUCCEEDED;
  };
  @action
  onAggCaseDataErr = () => {
    this.aggCaseDataState = ASYNC_STATES.FAILED;
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
    this.activeTab = tab;
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

export default ObservableUIStore;

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

class ObservableUIStore {
  @observable sidebarOpen = false;
  @observable sidebarSelectedGroupKeys = [];
  @observable caseDataState = ASYNC_STATES.STARTED;
  @observable aggCaseDataState = ASYNC_STATES.STARTED;
  // @observable activeTab = TABS.TAB_GROUP;
  @observable activeTab = TABS.TAB_EXAMPLE;

  @observable keysPressed = [];

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

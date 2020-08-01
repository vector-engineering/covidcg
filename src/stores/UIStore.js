import { observable, action, toJS } from 'mobx';

const STARTED = 'STARTED';
const SUCCEEDED = 'SUCCEEDED';
const FAILED = 'FAILED';

export const asyncStates = {
  STARTED,
  SUCCEEDED,
  FAILED,
};

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
  @observable caseDataState = STARTED;
  @observable aggCaseDataState = STARTED;
  @observable activeTab = 'group';

  @observable keysPressed = [];

  @action
  onCaseDataStateStarted = () => {
    this.caseDataState = STARTED;
  };

  @action
  onCaseDataStateFinished = () => {
    this.caseDataState = SUCCEEDED;
  };

  @action
  onCaseDataStateErr = () => {
    this.caseDataState = FAILED;
  };

  @action
  onAggCaseDataStarted = () => {
    this.aggCaseDataState = STARTED;
  };

  @action
  onAggCaseDataFinished = () => {
    this.aggCaseDataState = SUCCEEDED;
  };

  @action
  onAggCaseDataErr = () => {
    this.aggCaseDataState = FAILED;
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

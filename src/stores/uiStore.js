import { observable, action } from 'mobx';

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

class UiStore {
  @observable sidebarOpen = false;
  @observable sidebarSelectedGroupKeys = [];
  @observable dataState = STARTED;

  @action
  onDataChangeStarted = () => {
    this.dataState = STARTED;
  };

  @action
  onDataChangeFinished = () => {
    this.dataState = SUCCEEDED;
  };

  @action
  onDataChangeErr = () => {
    this.dataState = FAILED;
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
}

export default UiStore;

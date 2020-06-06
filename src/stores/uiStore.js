import { observable, action } from 'mobx';

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

  @action
  openSidebar = () => {
    console.log('yo');
    this.sidebarOpen = true;
  };

  @action
  closeSidebar = () => {
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

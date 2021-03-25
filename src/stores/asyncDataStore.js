import { hostname } from '../config';
import { observable, action, runInAction } from 'mobx';
import { ASYNC_STATES } from '../constants/defs.json';

class ObservableAsyncDataStore {
  @observable status = ASYNC_STATES.UNINITIALIZED;
  data = {};

  @action
  async fetchData() {
    this.status = ASYNC_STATES.STARTED;
    fetch(hostname + '/init')
      .then((res) => {
        if (!res.ok) {
          throw res;
        }
        return res.json();
      })
      .then((data) => {
        runInAction(() => {
          this.data = data;
          this.status = ASYNC_STATES.SUCCEEDED;
        });
      })
      .catch((err) => {
        runInAction(() => {
          this.status = ASYNC_STATES.FAILED;
        });
        err.text().then((errMsg) => {
          console.error('Error getting initial data');
          console.error(errMsg);
        });
      });
  }
}

export default ObservableAsyncDataStore;

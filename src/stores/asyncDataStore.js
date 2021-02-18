import { init_endpoint } from '../config';
import { observable, action, runInAction } from 'mobx';
import { ASYNC_STATES } from '../constants/defs.json';

class ObservableAsyncDataStore {
  @observable status = ASYNC_STATES.UNINITIALIZED;
  data = {};

  @action
  async fetchData() {
    this.status = ASYNC_STATES.STARTED;
    try {
      const res = await fetch(init_endpoint);
      const data = await res.json();
      runInAction(() => {
        this.data = data;
        this.status = ASYNC_STATES.SUCCEEDED;
      });
    } catch (e) {
      console.error(e);
      runInAction(() => {
        this.status = ASYNC_STATES.FAILED;
      });
    }
  }
}

export default ObservableAsyncDataStore;

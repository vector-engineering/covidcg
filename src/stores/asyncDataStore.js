import { observable, action, toJS } from 'mobx';
import { ASYNC_STATES } from '../constants/UI';

class ObservableAsyncDataStore {
  @observable status = ASYNC_STATES.UNINITIALIZED;
  @observable data = {};

  @action
  requestData() {
    this.status = ASYNC_STATES.STARTED;
    // setTimeout(() => {
    //   this.status = ASYNC_STATES.SUCCEEDED;
    // }, 2000);
  }
}

export default ObservableAsyncDataStore;

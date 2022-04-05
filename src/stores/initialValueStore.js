import { observable, action, runInAction } from 'mobx';

import { ASYNC_STATES } from '../constants/defs.json';

class InitialValueStore {
  @observable configStore = {};
  @observable plotSettingsStore = {};
  @observable groupDataStore = {};

  @observable status = ASYNC_STATES.UNINITIALIZED;

  constructor() {
    this.getInitialValues();
  }

  @action
  async getInitialValues() {
    this.status = ASYNC_STATES.STARTED;
    import('../config')
      .then((configMod) => {
        return configMod.config;
      })
      .then((config) => {
        import(`../constants/initialValues.${config.virus}.js`)
          .then((mod) => {
            return mod.default();
          })
          .then((allInitialValues) => {
            runInAction(() => {
              this.configStore = allInitialValues.configStore;
              this.plotSettingsStore = allInitialValues.plotSettingsStore;
              this.groupDataStore = allInitialValues.groupDataStore;
              this.status = ASYNC_STATES.SUCCEEDED;
            });
          });
      })
      .catch((err) => {
        runInAction(() => {
          this.status = ASYNC_STATES.FAILED;
        });
        let prefix = 'Error setting initial values';
        if (!(typeof err.text === 'function')) {
          console.error(prefix, err);
        } else {
          err.text().then((errMsg) => {
            console.error(prefix, errMsg);
          });
        }
      });
  }
}

export default InitialValueStore;

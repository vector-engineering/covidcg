import { observable, action, runInAction } from 'mobx';

import { ASYNC_STATES } from '../constants/defs.json';

export default class InitialValueStore {
  @observable configStore = {};
  @observable plotSettingsStore = {};

  @observable status = ASYNC_STATES.UNINITIALIZED;

  constructor() {
    this.getInitialValues();
  }

  init() {}

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
            console.log(mod);
            return mod.default();
          })
          .then((allInitialValues) => {
            runInAction(() => {
              this.configStore = allInitialValues.configStore;
              this.plotSettingsStore = allInitialValues.plotSettingsStore;
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

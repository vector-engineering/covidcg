import { observable, action, runInAction } from 'mobx';

export class InitialValueStore {
  @observable configStore = {};

  constructor() {}

  init() {
    this.getInitialValues();
  }

  @action
  async getInitialValues() {
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
            });
          });
      });
  }
}

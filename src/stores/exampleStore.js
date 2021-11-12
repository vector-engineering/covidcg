import { observable, action, runInAction } from 'mobx';
import { rootStoreInstance } from './rootStore';

export class ExampleStore {
  @observable examples = [];

  constructor() {}

  init() {
    this.getExamples();
  }

  @action
  async getExamples() {
    import('../config')
      .then((configMod) => {
        return configMod.config;
      })
      .then((config) => {
        import(`../components/Example/examples.${config.virus}.js`)
          .then((mod) => {
            return mod.default(rootStoreInstance.locationDataStore.selectTree);
          })
          .then((examplesList) => {
            runInAction(() => {
              this.examples = examplesList;
            });
          });
      });
  }
}

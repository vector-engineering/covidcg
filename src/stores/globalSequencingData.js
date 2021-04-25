import { action } from 'mobx';
import { hostname } from '../config';
import { rootStoreInstance } from './rootStore';
import { TABS } from '../constants/defs.json';

export class GlobalSequencingDataStore {
  countryScoreData;

  constructor() {}

  init() {
    this.countryScoreData = [];

    if (rootStoreInstance.UIStore.activeTab === TABS.TAB_GLOBAL_SEQUENCES) {
      this.fetchGlobalSequencingData();
    }
  }

  @action
  async fetchGlobalSequencingData() {
    rootStoreInstance.UIStore.onGlobalSequencingDataStarted();
    fetch(hostname + '/country_score', {
      method: 'GET',
      headers: {
        Accept: 'application/json',
      },
    })
      .then((res) => {
        if (!res.ok) {
          throw res;
        }
        return res.json();
      })
      .then((pkg) => {
        this.countryScoreData = pkg.country_score;
        rootStoreInstance.UIStore.onGlobalSequencingDataFinished();
      })
      .catch((err) => {
        err.text().then((errMsg) => {
          console.error(errMsg);
          rootStoreInstance.UIStore.onGlobalSequencingDataErr();
        });
      });
  }
}

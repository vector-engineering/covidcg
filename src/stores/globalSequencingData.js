// import { observable, toJS, action } from 'mobx';
import { asyncDataStoreInstance } from '../components/App';

export class GlobalSequencingDataStore {
  countryScoreData;

  constructor() {}

  init() {
    this.countryScoreData = asyncDataStoreInstance.data.country_score;
  }
}

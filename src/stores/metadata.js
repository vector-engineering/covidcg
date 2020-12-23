// import { observable, toJS, action } from 'mobx';
import { asyncDataStoreInstance } from '../components/App';

export class MetadataStore {
  metadataMap;

  constructor() {}

  init() {
    this.metadataMap = asyncDataStoreInstance.data.metadata_map;
  }

  getMetadataValueFromId(field, id) {
    return this.metadataMap[field][id];
  }
}

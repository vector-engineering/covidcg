import { action, toJS } from 'mobx';
import { hostname } from '../config';
import { asyncDataStoreInstance } from '../components/App';
import { rootStoreInstance } from './rootStore';
import { ASYNC_STATES } from '../constants/defs.json';

export class MetadataStore {
  UIStoreInstance;
  dataStoreInstance;

  metadataMap;

  constructor() {}

  init() {
    this.UIStoreInstance = rootStoreInstance.UIStore;
    this.dataStoreInstance = rootStoreInstance.dataStore;

    this.metadataMap = asyncDataStoreInstance.data.metadata_map;
  }

  @action
  async fetchMetadataFields() {
    this.UIStoreInstance.onMetadataFieldStarted();
    const metadataCounts = toJS(this.dataStoreInstance.metadataCounts);
    const requestedKeys = {};

    Object.keys(metadataCounts).forEach((field) => {
      requestedKeys[field] = Object.keys(metadataCounts[field]).map((key) =>
        parseInt(key)
      );
    });

    const res = await fetch(hostname + '/metadata_fields', {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(requestedKeys),
    });
    this.metadataMap = await res.json();

    this.UIStoreInstance.onMetadataFieldFinished();
  }

  getMetadataValueFromId(field, id) {
    if (this.UIStoreInstance.metadataFieldState !== ASYNC_STATES.SUCCEEDED) {
      return 'Waiting...';
    }

    return this.metadataMap[field][id];
  }
}

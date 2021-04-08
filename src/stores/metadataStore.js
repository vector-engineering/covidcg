import { action, toJS } from 'mobx';
import { config, hostname } from '../config';
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

    this.metadataMap = {};
    Object.keys(config.metadata_cols).forEach((field) => {
      this.metadataMap[field] = {};
    });
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

    if (Object.keys(requestedKeys).length === 0) {
      this.UIStoreInstance.onMetadataFieldFinished();
      return;
    }

    fetch(hostname + '/metadata_fields', {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(requestedKeys),
    })
      .then((res) => {
        if (!res.ok) {
          throw res;
        }
        return res.json();
      })
      .then((metadataMap) => {
        this.metadataMap = metadataMap;
        this.UIStoreInstance.onMetadataFieldFinished();
      })
      .catch((err) => {
        err.text().then((errMsg) => {
          console.error(errMsg);
          this.UIStoreInstance.onMetadataFieldErr();
        });
      });
  }

  getMetadataValueFromId(field, id) {
    if (this.UIStoreInstance.metadataFieldState === ASYNC_STATES.SUCCEEDED) {
      if (!Object.prototype.hasOwnProperty.call(this.metadataMap, field)) {
        return undefined;
      }

      return this.metadataMap[field][id];
    } else if (
      this.UIStoreInstance.metadataFieldState === ASYNC_STATES.STARTED ||
      this.UIStoreInstance.metadataFieldState === ASYNC_STATES.UNINITIALIZED
    ) {
      return 'Waiting...';
    } else {
      return 'Error fetching fields';
    }
  }
}

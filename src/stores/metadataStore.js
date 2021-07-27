import { action, toJS } from 'mobx';
import { config, hostname } from '../config';
import { rootStoreInstance } from './rootStore';
import { ASYNC_STATES } from '../constants/defs.json';

import { getLocationIdsByNode } from '../utils/location';

export class MetadataStore {
  // Map of metadata value ID --> metadata value string
  metadataMap;
  // Map of metadata value ID --> counts
  metadataCounts;

  constructor() {}

  init() {
    this.metadataMap = new Map();
    this.metadataCounts = new Map();

    Object.keys(config.metadata_cols).forEach((field) => {
      this.metadataMap.set(field, new Map());
      this.metadataCounts.set(field, new Map());
    });
  }

  @action
  async fetchMetadataFields() {
    rootStoreInstance.UIStore.onMetadataFieldStarted();

    fetch(hostname + '/metadata', {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        start_date: toJS(rootStoreInstance.configStore.startDate),
        end_date: toJS(rootStoreInstance.configStore.endDate),
        location_ids: getLocationIdsByNode(
          toJS(rootStoreInstance.configStore.selectedLocationNodes)
        ),
        selected_metadata_fields:
          rootStoreInstance.configStore.getSelectedMetadataFields(),
      }),
    })
      .then((res) => {
        if (!res.ok) {
          throw res;
        }
        return res.json();
      })
      .then((metadataRecords) => {
        // Each row in this array of records is structured as:
        // { "field", "val_id", "count", "val_str" }
        metadataRecords.forEach((record) => {
          this.metadataMap.get(record.field).set(record.val_id, record.val_str);
          this.metadataCounts
            .get(record.field)
            .set(record.val_id, record.count);
        });

        rootStoreInstance.UIStore.onMetadataFieldFinished();
      })
      .catch((err) => {
        if (!(typeof err.text === 'function')) {
          console.error(err);
        } else {
          err.text().then((errMsg) => {
            console.error(errMsg);
          });
        }
        rootStoreInstance.UIStore.onMetadataFieldErr();
      });
  }

  getMetadataValueFromId(field, id) {
    if (
      rootStoreInstance.UIStore.metadataFieldState === ASYNC_STATES.SUCCEEDED
    ) {
      if (!this.metadataMap.has(field)) {
        return undefined;
      }

      return this.metadataMap.get(field).get(id);
    } else if (
      rootStoreInstance.UIStore.metadataFieldState === ASYNC_STATES.STARTED ||
      rootStoreInstance.UIStore.metadataFieldState ===
        ASYNC_STATES.UNINITIALIZED
    ) {
      return 'Waiting...';
    } else {
      return 'Error fetching fields';
    }
  }

  getMetadataCountsFromId(field, id) {
    if (
      rootStoreInstance.UIStore.metadataFieldState === ASYNC_STATES.SUCCEEDED
    ) {
      if (!this.metadataCounts.has(field)) {
        return undefined;
      }

      return this.metadataCounts.get(field).get(id);
    } else if (
      rootStoreInstance.UIStore.metadataFieldState === ASYNC_STATES.STARTED ||
      rootStoreInstance.UIStore.metadataFieldState ===
        ASYNC_STATES.UNINITIALIZED
    ) {
      return null;
    } else {
      return null;
    }
  }
}

/**
 * Data for groups (i.e., lineages/clades)
 */

import { action } from 'mobx';
import { hostname } from '../config';
import { rootStoreInstance } from './rootStore';

import { downloadBlobURL } from '../utils/download';

export class GroupDataStore {
  groupSnvFrequency;

  constructor() {
    this.groupSnvFrequency = {};
  }

  init() {}

  @action
  async fetchGroupSnvFrequencyData({ group, snvType, consensusThreshold }) {
    rootStoreInstance.UIStore.onGroupSnvFrequencyStarted();
    return fetch(hostname + '/group_snv_frequencies', {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        group,
        snv_type: snvType,
        consensus_threshold: consensusThreshold,
      }),
    })
      .then((res) => {
        if (!res.ok) {
          throw res;
        }
        return res.json();
      })
      .then((pkg) => {
        if (
          !Object.prototype.hasOwnProperty.call(this.groupSnvFrequency, group)
        ) {
          this.groupSnvFrequency[group] = {};
        }
        this.groupSnvFrequency[group][snvType] = pkg;

        rootStoreInstance.UIStore.onGroupSnvFrequencyFinished();

        return pkg;
      })
      .catch((err) => {
        err.text().then((errMsg) => {
          console.error(errMsg);
          rootStoreInstance.UIStore.onGroupSnvFrequencyErr();
        });
      });
  }

  async downloadGroupSnvFrequencyData({ group, snvType, consensusThreshold }) {
    rootStoreInstance.UIStore.onDownloadStarted();
    this.fetchGroupSnvFrequencyData({
      group,
      snvType,
      consensusThreshold,
    }).then((pkg) => {
      rootStoreInstance.UIStore.onDownloadFinished();
      const blob = new Blob([JSON.stringify(pkg)]);
      const url = URL.createObjectURL(blob);
      downloadBlobURL(url, 'consensus_mutations.json');
    });
  }
}

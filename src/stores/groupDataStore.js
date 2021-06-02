/**
 * Data for groups (i.e., lineages/clades)
 */

import { action, observable } from 'mobx';
import { config, hostname } from '../config';
import { rootStoreInstance } from './rootStore';
import { asyncDataStoreInstance } from '../components/App';

import { downloadBlobURL } from '../utils/download';

export class GroupDataStore {
  // Actively selected group type in the report tab
  @observable activeGroupType;
  @observable selectedGroups;
  @observable groupSnvType;
  @observable consensusThreshold;

  groups;
  @observable groupSelectTree;
  groupSnvFrequency;

  constructor() {
    this.activeGroupType = Object.keys(config['group_cols'])[0];
    this.selectedGroups = ['B.1.1.7', 'B.1.351', 'P.2'];
    this.groupSnvType = 'gene_aa';
    this.consensusThreshold = 0.7;

    this.groupSnvFrequency = {};

    this.groups = {};
    this.groupSelectTree = {};

    Object.keys(config['group_cols']).forEach((group) => {
      this.groups[group] = [];
      this.groupSelectTree[group] = [];
    });
  }

  init() {
    // Load all groups from the server
    asyncDataStoreInstance.data.groups.forEach((record) => {
      let groupName = record['group'];
      delete record['group'];
      this.groups[groupName].push(record);
    });

    // Construct selection trees
    Object.keys(this.groups).forEach((groupName) => {
      this.groups[groupName].forEach((group) => {
        this.groupSelectTree[groupName].push({
          label: group.name,
          value: group.name,
          checked: false,
        });
      });
    });
  }

  @action
  updateActiveGroupType(activeGroupType) {
    this.activeGroupType = activeGroupType;
  }

  @action
  updateGroupSnvType(groupSnvType) {
    this.groupSnvType = groupSnvType;
  }

  @action
  updateConsensusThreshold(consensusThreshold) {
    this.consensusThreshold = consensusThreshold;
  }

  @action
  updateSelectedGroups(selectedGroups) {
    this.selectedGroups = selectedGroups;
  }

  getActiveGroupTypePrettyName() {
    return config.group_cols[this.activeGroupType].title;
  }

  getGroupSnvTypePrettyName() {
    if (this.groupSnvType === 'dna') {
      return 'NT';
    } else if (this.groupSnvType === 'gene_aa') {
      return 'Gene AA';
    } else if (this.groupSnvType === 'protein_aa') {
      return 'Protein AA';
    }
  }

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

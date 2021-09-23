/**
 * Data for groups (i.e., lineages/clades)
 */

import { action, observable } from 'mobx';
import { config, hostname } from '../config';
import { rootStoreInstance } from './rootStore';
import { asyncDataStoreInstance } from '../components/App';
import { GROUPS } from '../constants/defs.json';

import { downloadBlobURL } from '../utils/download';
import { mutationHeatmapToPymolScript } from '../utils/pymol';

export const initialValues = {
  activeGroupType: Object.keys(config['group_cols'])[0],
  selectedGroups: ['AY.4', 'B.1.617.2', 'B.1.1.7', 'B.1.351', 'P.2'],
  groupSnvType: 'protein_aa',
};

export class GroupDataStore {
  // Actively selected group type in the report tab
  @observable activeGroupType = initialValues.activeGroupType;
  @observable selectedGroups = initialValues.selectedGroups;
  @observable groupSnvType = initialValues.groupSnvType;

  groups;
  @observable groupSelectTree;
  groupSnvFrequency;

  constructor() {
    this.groupSnvFrequency = {};
    // Provided by the server
    // Array of records { name: string, color: string }
    this.groups = {};
    this.groupSelectTree = {};
    this.groupColors = {};

    Object.keys(config['group_cols']).forEach((group) => {
      this.groups[group] = [];
      this.groupSelectTree[group] = [];
      this.groupColors[group] = {};
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
          checked: this.selectedGroups.includes(group.name),
        });
      });
    });

    // Assign group colors
    Object.keys(this.groups).forEach((groupName) => {
      // Assign color for the "Other" group
      this.groupColors[groupName][GROUPS.OTHER_GROUP] = '#AAA';
      this.groups[groupName].forEach((group) => {
        this.groupColors[groupName][group.name] = group.color;
      });
    });
  }

  hasGroupFrequencyData() {
    if (
      !Object.prototype.hasOwnProperty.call(
        this.groupSnvFrequency,
        this.activeGroupType
      )
    ) {
      return false;
    } else if (
      !Object.prototype.hasOwnProperty.call(
        this.groupSnvFrequency[this.activeGroupType],
        this.groupSnvType
      )
    ) {
      return false;
    } else {
      return true;
    }
  }

  @action
  updateActiveGroupType = (activeGroupType) => {
    this.activeGroupType = activeGroupType;

    // If we don't have the data for this combo yet, then fetch it now
    if (!this.hasGroupFrequencyData()) {
      this.fetchGroupSnvFrequencyData({
        group: this.activeGroupType,
        snvType: this.groupSnvType,
        consensusThreshold: 0,
      });
    }
  };

  @action
  updateGroupSnvType = (groupSnvType) => {
    this.groupSnvType = groupSnvType;

    // If we don't have the data for this combo yet, then fetch it now
    if (!this.hasGroupFrequencyData()) {
      this.fetchGroupSnvFrequencyData({
        group: this.activeGroupType,
        snvType: this.groupSnvType,
        consensusThreshold: 0,
      });
    }

    // Hide ORF1a in NT/gene_aa mode by default
    if (groupSnvType === 'dna' || groupSnvType === 'gene_aa') {
      rootStoreInstance.plotSettingsStore.setReportMutationListHidden([
        'ORF1a',
      ]);
    }
    // Otherwise clear the hidden list
    else {
      rootStoreInstance.plotSettingsStore.setReportMutationListHidden([]);
    }
  };

  @action
  updateSelectedGroups = (selectedGroups) => {
    this.selectedGroups = selectedGroups;

    // Update the groupSelectTree as well
    const selectTree = JSON.parse(JSON.stringify(this.groupSelectTree));
    selectTree[this.activeGroupType].forEach((group) => {
      if (this.selectedGroups.includes(group.value)) {
        group.checked = true;
      } else {
        group.checked = false;
      }
    });
    this.groupSelectTree = selectTree;

    // Update the selected active group from the structural viewer
    // if we removed the current structural active group
    if (
      !this.selectedGroups.includes(
        rootStoreInstance.plotSettingsStore.reportStructureActiveGroup
      )
    ) {
      // Set it to the first selected group
      rootStoreInstance.plotSettingsStore.setReportStructureActiveGroup(
        this.selectedGroups[0]
      );
    }
  };

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

  getGroupColor(groupKey, group) {
    return this.groupColors[groupKey][group];
  }

  @action
  async fetchGroupSnvFrequencyData({ group, snvType, consensusThreshold }) {
    // Skip the download if we already have the requested data
    if (
      Object.prototype.hasOwnProperty.call(this.groupSnvFrequency, group) &&
      Object.prototype.hasOwnProperty.call(
        this.groupSnvFrequency[group],
        snvType
      )
    ) {
      // eslint-disable-next-line no-unused-vars
      return new Promise((resolve, reject) => {
        resolve(this.groupSnvFrequency[group][snvType]);
      });
    }

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
        if (!(typeof err.text === 'function')) {
          console.error(err);
        } else {
          err.text().then((errMsg) => {
            console.error(errMsg);
          });
        }
        rootStoreInstance.UIStore.onGroupSnvFrequencyErr();
      });
  }

  @action
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

  getStructureSnvs() {
    return this.groupSnvFrequency[this.activeGroupType]['protein_aa'].filter(
      (groupSnv) =>
        groupSnv.name ===
          rootStoreInstance.plotSettingsStore.reportStructureActiveGroup &&
        groupSnv.protein ===
          rootStoreInstance.plotSettingsStore.reportStructureActiveProtein
    );
  }

  downloadStructureMutationData() {
    const blob = new Blob([JSON.stringify(this.getStructureSnvs())]);
    const url = URL.createObjectURL(blob);
    downloadBlobURL(
      url,
      `mutations_${rootStoreInstance.plotSettingsStore.reportStructureActiveProtein}_${rootStoreInstance.plotSettingsStore.reportStructureActiveGroup}.json`
    );
  }

  downloadStructurePymolScript(opts) {
    const script = mutationHeatmapToPymolScript({
      activeProtein:
        rootStoreInstance.plotSettingsStore.reportStructureActiveProtein,
      activeGroup:
        rootStoreInstance.plotSettingsStore.reportStructureActiveGroup,
      pdbId: rootStoreInstance.plotSettingsStore.reportStructurePdbId,
      snvs: this.getStructureSnvs(),
      ...opts,
    });
    const blob = new Blob([script]);
    const url = URL.createObjectURL(blob);
    downloadBlobURL(
      url,
      `heatmap_${rootStoreInstance.plotSettingsStore.reportStructureActiveProtein}_${rootStoreInstance.plotSettingsStore.reportStructureActiveGroup}.py`
    );
  }
}

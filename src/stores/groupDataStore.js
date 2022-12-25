/**
 * Data for groups (i.e., lineages/clades)
 */

import { action, observable } from 'mobx';
import { config, hostname } from '../config';
import { rootStoreInstance } from './rootStore';
import { asyncDataStoreInstance } from '../components/App';
import { groupDataStore as initialGroupDataStore } from '../constants/initialValues';
import { GROUPS, PYMOL_SCRIPT_TYPES } from '../constants/defs.json';

import { downloadBlobURL } from '../utils/download';

import { getProtein } from '../utils/gene_protein';
import { savePymolScript } from '../utils/pymol';
import { updateURLFromParams } from '../utils/updateQueryParam';

export class GroupDataStore {
  initialValues = {};
  // Actively selected group type in the report tab
  @observable activeGroupType = '';
  @observable selectedReportGroups = [];
  @observable groupMutationType = '';
  @observable activeReference = '';

  groups;
  @observable groupSelectTree;
  groupMutationFrequency;

  constructor() {
    this.groupMutationFrequency = {};
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
    // Load intial values
    this.initialValues = initialGroupDataStore;
    Object.keys(this.initialValues).forEach((key) => {
      this[key] = this.initialValues[key];
    });
    // Load all groups from the server
    asyncDataStoreInstance.data.groups.forEach((record) => {
      let groupName = record['group'];
      //delete record['group'];
      this.groups[groupName].push(record);
    });

    // Construct selection trees
    Object.keys(this.groups).forEach((groupName) => {
      this.groups[groupName].forEach((group) => {
        this.groupSelectTree[groupName].push({
          label: group.name,
          value: group.name,
          checked: this.selectedReportGroups.includes(group.name),
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

  @action
  applyPendingChanges = (pending, fetch = true) => {
    // Overwrite any of our fields here with the pending ones
    Object.keys(pending).forEach((field) => {
      this[field] = pending[field];
    });

    // Hide ORF1a in NT/gene_aa mode by default
    if (
      pending.groupMutationType === 'dna' ||
      pending.groupMutationType === 'gene_aa'
    ) {
      rootStoreInstance.plotSettingsStore.setReportMutationListHidden([
        'ORF1a',
      ]);
    }
    // Otherwise clear the hidden list
    else {
      rootStoreInstance.plotSettingsStore.setReportMutationListHidden([]);
    }

    // Update the groupSelectTree as well
    const selectTree = JSON.parse(JSON.stringify(this.groupSelectTree));
    selectTree[this.activeGroupType].forEach((group) => {
      if (this.selectedReportGroups.includes(group.value)) {
        group.checked = true;
      } else {
        group.checked = false;
      }
    });
    this.groupSelectTree = selectTree;

    // Update the selected active group from the structural viewer
    // if we removed the current structural active group
    if (
      !this.selectedReportGroups.includes(
        rootStoreInstance.plotSettingsStore.reportStructureActiveGroup
      )
    ) {
      // Set it to the first selected group
      rootStoreInstance.plotSettingsStore.setReportStructureActiveGroup(
        this.selectedReportGroups[0]
      );
    }

    // If we don't have the data for this combo yet, then fetch it now
    if (fetch || !this.hasGroupFrequencyData()) {
      this.fetchGroupMutationFrequencyData({
        group: this.activeGroupType,
        mutationType: this.groupMutationType,
        consensusThreshold: 0,
        selectedReference: this.activeReference,
      });
    }

    this.updateURL(pending);
  };

  updateURL = (pending) => {
    const urlParams = rootStoreInstance.urlMonitor.urlParams;
    Object.keys(pending).forEach((field) => {
      if (field === 'selectedReportGroups') {
        urlParams.set(field, pending[field].join(','));
      } else {
        urlParams.set(field, String(pending[field]));
      }

      if (pending[field] === this.initialValues[field]) {
        // Only display non-default fields in the url
        urlParams.delete(field);
      }
    });

    // Update URL
    updateURLFromParams(urlParams);

    rootStoreInstance.urlMonitor.urlParams = urlParams;
  };

  hasGroupFrequencyData() {
    if (
      !Object.prototype.hasOwnProperty.call(
        this.groupMutationFrequency,
        this.activeGroupType
      )
    ) {
      return false;
    } else if (
      !Object.prototype.hasOwnProperty.call(
        this.groupMutationFrequency[this.activeGroupType],
        this.groupMutationType
      )
    ) {
      return false;
    } else if (
      !Object.prototype.hasOwnProperty.call(
        this.groupMutationFrequency[this.activeGroupType][
          this.groupMutationType
        ],
        '0'
      )
    ) {
      return false;
    } else {
      return true;
    }
  }

  getActiveGroupTypePrettyName() {
    return config.group_cols[this.activeGroupType].title;
  }

  getGroupMutationTypePrettyName() {
    if (this.groupMutationType === 'dna') {
      return 'NT';
    } else if (this.groupMutationType === 'gene_aa') {
      return 'Gene AA';
    } else if (this.groupMutationType === 'protein_aa') {
      return 'Protein AA';
    }
  }

  getGroupColor(groupKey, group) {
    return this.groupColors[groupKey][group];
  }

  @action
  async fetchGroupMutationFrequencyData({
    group,
    mutationType,
    consensusThreshold,
    selectedReference,
  }) {
    // Skip the download if we already have the requested data
    if (
      Object.prototype.hasOwnProperty.call(
        this.groupMutationFrequency,
        group
      ) &&
      Object.prototype.hasOwnProperty.call(
        this.groupMutationFrequency[group],
        mutationType
      ) &&
      Object.prototype.hasOwnProperty.call(
        this.groupMutationFrequency[group][mutationType],
        consensusThreshold.toString()
      )
    ) {
      // eslint-disable-next-line no-unused-vars
      return new Promise((resolve, reject) => {
        resolve(
          this.groupMutationFrequency[group][mutationType][
            consensusThreshold.toString()
          ]
        );
      });
    }

    rootStoreInstance.UIStore.onGroupMutationFrequencyStarted();
    return fetch(hostname + '/group_mutation_frequencies', {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        group_key: group,
        mutation_type: mutationType,
        consensus_threshold: consensusThreshold,
        selected_reference: selectedReference,
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
          !Object.prototype.hasOwnProperty.call(
            this.groupMutationFrequency,
            group
          )
        ) {
          this.groupMutationFrequency[group] = {};
        }
        if (
          !Object.prototype.hasOwnProperty.call(
            this.groupMutationFrequency[group],
            mutationType
          )
        ) {
          this.groupMutationFrequency[group][mutationType] = {};
        }
        this.groupMutationFrequency[group][mutationType][
          consensusThreshold.toString()
        ] = pkg;

        rootStoreInstance.UIStore.onGroupMutationFrequencyFinished();

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
        rootStoreInstance.UIStore.onGroupMutationFrequencyErr();
      });
  }

  @action
  async downloadGroupMutationFrequencyData({
    group,
    mutationType,
    consensusThreshold,
    selectedReference,
  }) {
    rootStoreInstance.UIStore.onDownloadStarted();
    this.fetchGroupMutationFrequencyData({
      group,
      mutationType,
      consensusThreshold,
      selectedReference,
    }).then((pkg) => {
      rootStoreInstance.UIStore.onDownloadFinished();
      const blob = new Blob([JSON.stringify(pkg)]);
      const url = URL.createObjectURL(blob);
      downloadBlobURL(url, 'consensus_mutations.json');
    });
  }

  getStructureMutations() {
    return this.groupMutationFrequency[this.activeGroupType]['protein_aa'][
      '0'
    ].filter(
      (groupMutation) =>
        groupMutation.name ===
          rootStoreInstance.plotSettingsStore.reportStructureActiveGroup &&
        groupMutation.feature ===
          rootStoreInstance.plotSettingsStore.reportStructureActiveProtein
    );
  }

  downloadStructureMutationData() {
    const blob = new Blob([JSON.stringify(this.getStructureMutations())]);
    const url = URL.createObjectURL(blob);
    downloadBlobURL(
      url,
      `mutations_${rootStoreInstance.plotSettingsStore.reportStructureActiveProtein}_${rootStoreInstance.plotSettingsStore.reportStructureActiveGroup}.json`
    );
  }

  downloadStructurePymolScript(opts) {
    let filename;
    if (opts.scriptType === PYMOL_SCRIPT_TYPES.COMMANDS) {
      filename = `heatmap_${rootStoreInstance.plotSettingsStore.reportStructureActiveProtein}_${rootStoreInstance.plotSettingsStore.reportStructureActiveGroup}.txt`;
    } else if (opts.scriptType === PYMOL_SCRIPT_TYPES.SCRIPT) {
      filename = `heatmap_${rootStoreInstance.plotSettingsStore.reportStructureActiveProtein}_${rootStoreInstance.plotSettingsStore.reportStructureActiveGroup}.py`;
    }

    savePymolScript({
      opts,
      filename,
      activeProtein: getProtein(
        rootStoreInstance.plotSettingsStore.reportStructureActiveProtein,
        rootStoreInstance.configStore.selectedReference // TODO: track reference for report in plotSettingsStore?
      ),
      pdbId: rootStoreInstance.plotSettingsStore.reportStructurePdbId,
      proteinStyle:
        rootStoreInstance.plotSettingsStore.reportStructureProteinStyle,
      assemblies: rootStoreInstance.plotSettingsStore.reportStructureAssemblies,
      activeAssembly:
        rootStoreInstance.plotSettingsStore.reportStructureActiveAssembly,
      entities: rootStoreInstance.plotSettingsStore.reportStructureEntities,
      mutations: this.getStructureMutations(),
      mutationColorField: 'fraction',
    });
  }
}

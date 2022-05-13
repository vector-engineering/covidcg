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
import {
  mutationHeatmapToPymolScript,
  mutationHeatmapToPymolCommands,
} from '../utils/pymol';

export class GroupDataStore {
  initialValues = {};
  // Actively selected group type in the report tab
  @observable activeGroupType = '';
  @observable selectedGroups = [];
  @observable groupMutationType = '';

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

  @action
  updateActiveGroupType = (activeGroupType) => {
    this.activeGroupType = activeGroupType;

    // If we don't have the data for this combo yet, then fetch it now
    if (!this.hasGroupFrequencyData()) {
      this.fetchGroupMutationFrequencyData({
        group: this.activeGroupType,
        mutationType: this.groupMutationType,
        consensusThreshold: 0,
      });
    }
  };

  @action
  updateGroupMutationType = (groupMutationType) => {
    this.groupMutationType = groupMutationType;

    // If we don't have the data for this combo yet, then fetch it now
    if (!this.hasGroupFrequencyData()) {
      this.fetchGroupMutationFrequencyData({
        group: this.activeGroupType,
        mutationType: this.groupMutationType,
        consensusThreshold: 0,
      });
    }

    // Hide ORF1a in NT/gene_aa mode by default
    if (groupMutationType === 'dna' || groupMutationType === 'gene_aa') {
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
        group,
        mutation_type: mutationType,
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
  }) {
    rootStoreInstance.UIStore.onDownloadStarted();
    this.fetchGroupMutationFrequencyData({
      group,
      mutationType,
      consensusThreshold,
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
        groupMutation.protein ===
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
    let script, outfile;
    if (opts.scriptType === PYMOL_SCRIPT_TYPES.COMMANDS) {
      script = mutationHeatmapToPymolCommands({
        activeProtein:
          rootStoreInstance.plotSettingsStore.reportStructureActiveProtein,
        pdbId: rootStoreInstance.plotSettingsStore.reportStructurePdbId,
        proteinStyle:
          rootStoreInstance.plotSettingsStore.reportStructureProteinStyle,
        assemblies:
          rootStoreInstance.plotSettingsStore.reportStructureAssemblies,
        activeAssembly:
          rootStoreInstance.plotSettingsStore.reportStructureActiveAssembly,
        entities: rootStoreInstance.plotSettingsStore.reportStructureEntities,
        mutations: this.getStructureMutations(),
        ...opts,
      });
      outfile = `heatmap_${rootStoreInstance.plotSettingsStore.reportStructureActiveProtein}_${rootStoreInstance.plotSettingsStore.reportStructureActiveGroup}.txt`;
    } else if (opts.scriptType === PYMOL_SCRIPT_TYPES.SCRIPT) {
      script = mutationHeatmapToPymolScript({
        activeProtein:
          rootStoreInstance.plotSettingsStore.reportStructureActiveProtein,
        activeGroup:
          rootStoreInstance.plotSettingsStore.reportStructureActiveGroup,
        pdbId: rootStoreInstance.plotSettingsStore.reportStructurePdbId,
        proteinStyle:
          rootStoreInstance.plotSettingsStore.reportStructureProteinStyle,
        assemblies:
          rootStoreInstance.plotSettingsStore.reportStructureAssemblies,
        activeAssembly:
          rootStoreInstance.plotSettingsStore.reportStructureActiveAssembly,
        entities: rootStoreInstance.plotSettingsStore.reportStructureEntities,
        mutations: this.getStructureMutations(),
        ...opts,
      });
      outfile = `heatmap_${rootStoreInstance.plotSettingsStore.reportStructureActiveProtein}_${rootStoreInstance.plotSettingsStore.reportStructureActiveGroup}.py`;
    }
    const blob = new Blob([script]);
    const url = URL.createObjectURL(blob);
    downloadBlobURL(url, outfile);
  }
}

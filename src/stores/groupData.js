import refSeq from '../../static_data/reference.json';

import { cladeColorArray } from '../constants/colors';
import { GROUPS } from '../constants/defs.json';
import { asyncDataStoreInstance } from '../components/App';

export class GroupDataStore {
  // Group -> SNV data as maps
  groupSnvMap;

  // Group -> color map
  groupColorMap;

  constructor() {
    this.groupSnvMap = {};
    this.groupColorMap = {};
  }

  init() {
    this.groupSnvMap = asyncDataStoreInstance.data.group_consensus_snps;

    // Build colors
    Object.keys(this.groupSnvMap).forEach((groupDef) => {
      // Initialize color map, with the "Other" group as grey by default
      const colorMap = {};
      colorMap[GROUPS.OTHER_GROUP] = '#AAA';

      Object.keys(this.groupSnvMap[groupDef]).forEach((group, ind) => {
        colorMap[group] = cladeColorArray[ind % cladeColorArray.length];
      });

      // Store in the master list
      this.groupColorMap[groupDef] = colorMap;
    });
  }
  getReferenceSequence() {
    return refSeq['ref_seq'];
  }

  getGroup(groupKey, group) {
    return this.groupSnvMap[groupKey][group];
  }
  getGroupColor(groupKey, group) {
    return this.groupColorMap[groupKey][group];
  }
}

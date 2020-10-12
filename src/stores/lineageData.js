// import lineageSnpData from '../../data/lineage_snp.json';
// import cladeSnpData from '../../data/clade_snp.json';
import refSeq from '../../static_data/reference.json';

import _ from 'underscore';

import { warmColors, coolColors, cladeColorArray } from '../constants/colors';
import { GROUP_KEYS } from '../constants/config';
import { GROUPS } from '../constants/groups';
import { asyncDataStoreInstance } from '../components/App';

export class LineageDataStore {
  // Internal counters for generating lineage colors
  coolColorInd;
  warmColorInd;

  // Internal counter for generating clade colors
  cladeColorInd;

  // Group -> SNV data as arrays
  lineageSnpData;
  cladeSnpData;

  // Group -> SNV data as maps
  groupSnvMap;

  // Colormaps
  cladeColorMap;
  lineageColorMap;
  groupColorMap;

  constructor() {
    // Init internal counters
    this.coolColorInd = 0;
    this.warmColorInd = 0;
    this.cladeColorInd = 0;

    this.groupSnvMap = {};
    this.groupSnvMap['lineage'] = {};
    this.groupSnvMap['clade'] = {};

    this.lineageColorMap = {};
    this.cladeColorMap = {};
    this.groupColorMap = {
      lineage: this.lineageColorMap,
      clade: this.cladeColorMap,
    };
  }

  init() {
    this.lineageSnpData = asyncDataStoreInstance.data.lineage_snp;
    this.cladeSnpData = asyncDataStoreInstance.data.clade_snp;

    // TODO: Do this in python instead
    this.lineageSnpData.forEach((lineageObj) => {
      this.groupSnvMap['lineage'][lineageObj.lineage] = lineageObj;
    });
    this.cladeSnpData.forEach((cladeObj) => {
      this.groupSnvMap['clade'][cladeObj.clade] = cladeObj;
    });

    this.lineageColorMap[GROUPS.OTHER_GROUP] = '#AAA';
    this.lineageSnpData.forEach((lineageObj) => {
      this.lineageColorMap[lineageObj.lineage] = this._getLineageColor(
        lineageObj.lineage
      );
    });

    this.cladeColorMap[GROUPS.OTHER_GROUP] = '#AAA';
    this.cladeSnpData.forEach((cladeObj) => {
      this.cladeColorMap[cladeObj.clade] = this._getCladeColor(cladeObj.clade);
    });
  }

  _getLineageColor = _.memoize((group) => {
    let color;
    if (group.charAt(0) === 'A') {
      color = warmColors[this.warmColorInd++];

      if (this.warmColorInd === warmColors.length) {
        this.warmColorInd = 0;
      }
    } else if (group.charAt(0) === 'B') {
      color = coolColors[this.coolColorInd++];

      if (this.coolColorInd === coolColors.length) {
        this.coolColorInd = 0;
      }
    }
    return color;
  });

  _getCladeColor = _.memoize(() => {
    const color = cladeColorArray[this.cladeColorInd++];
    // If we're at the end, then loop back to the beginning
    if (this.cladeColorInd === cladeColorArray.length) {
      this.cladeColorInd = 0;
    }

    return color;
  });

  loadLineageSnp() {
    return this.lineageSnpData;
  }
  loadCladeSnp() {
    return this.cladeSnpData;
  }

  getReferenceSequence() {
    return refSeq['ref_seq'];
  }

  getLineageColor = (lineage) => {
    return this.lineageColorMap[lineage];
  };

  getCladeColor = (clade) => {
    return this.cladeColorMap[clade];
  };

  getGroup(groupKey, group) {
    let snpData;
    let snpDataKey;
    if (groupKey === GROUP_KEYS.GROUP_LINEAGE) {
      snpData = this.lineageSnpData;
      snpDataKey = 'lineage';
    } else if (groupKey === GROUP_KEYS.GROUP_CLADE) {
      snpData = this.cladeSnpData;
      snpDataKey = 'clade';
    }
    let findObj = {};
    findObj[snpDataKey] = group;

    return _.findWhere(snpData, findObj);
  }
}

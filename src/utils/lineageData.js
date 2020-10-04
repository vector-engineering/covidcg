// import lineageSnpData from '../../data/lineage_snp.json';
// import cladeSnpData from '../../data/clade_snp.json';
import refSeq from '../../static_data/reference.json';

import _ from 'underscore';
import { intToDnaSnv, intToGeneAaSnv, intToProteinAaSnv } from './snpData';

import { warmColors, coolColors, cladeColorArray } from '../constants/colors';
import { GROUP_KEYS } from '../constants/config';
import { GROUPS } from '../constants/groups';
import { asyncDataStoreInstance } from '../components/App';

class LineageDataStore {
  coolColorInd;
  warmColorInd;
  cladeColorInd;

  lineageSnpData;
  cladeSnpData;
  cladeColorMap;
  lineageColorMap;

  init() {
    this.coolColorInd = 0;
    this.warmColorInd = 0;
    this.cladeColorInd = 0;

    this.lineageSnpData = asyncDataStoreInstance.data.lineage_snp;
    this.cladeSnpData = asyncDataStoreInstance.data.clade_snp;

    this.lineageColorMap = {};
    this.lineageColorMap[GROUPS.OTHER_GROUP] = '#AAA';
    this.lineageSnpData.forEach((lineageObj) => {
      this.lineageColorMap[lineageObj.lineage] = this._getLineageColor(
        lineageObj.lineage
      );
    });

    this.cladeColorMap = {};
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

  getDnaSnpsFromGroup(groupKey, group) {
    let groupObj = this.getGroup(groupKey, group);
    if (groupObj === undefined) {
      return [];
    }

    let snpIds = groupObj.dna_snp_ids;
    snpIds = _.reject(snpIds, (snpId) => snpId === '');
    return _.map(snpIds, (snpId) => intToDnaSnv(parseInt(snpId)));
  }

  getGeneAaSnpsFromGroup(groupKey, group) {
    let groupObj = this.getGroup(groupKey, group);
    if (groupObj === undefined) {
      return [];
    }
    let snpIds = groupObj.gene_aa_snp_ids;
    snpIds = _.reject(snpIds, (snpId) => snpId === '');
    return _.map(snpIds, (snpId) => intToGeneAaSnv(parseInt(snpId)));
  }

  getProteinAaSnpsFromGroup(groupKey, group) {
    let groupObj = this.getGroup(groupKey, group);
    if (groupObj === undefined) {
      return [];
    }
    let snpIds = groupObj.protein_aa_snp_ids;
    snpIds = _.reject(snpIds, (snpId) => snpId === '');
    return _.map(snpIds, (snpId) => intToProteinAaSnv(parseInt(snpId)));
  }
}

export const lineageDataStoreInstance = new LineageDataStore();

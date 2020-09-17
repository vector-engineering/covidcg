// import lineageSnpData from '../../data/lineage_snp.json';
// import cladeSnpData from '../../data/clade_snp.json';
import refSeq from '../../static_data/reference.json';

import _ from 'underscore';
import { intToDnaSnv, intToGeneAaSnv, intToProteinAaSnv } from './snpData';

import { warmColors, coolColors, cladeColorArray } from '../constants/colors';
import { GROUP_KEYS } from '../constants/config';
import { GROUPS } from '../constants/groups';

const lineageSnpData = {};
const cladeSnpData = {};

export function loadLineageSnp() {
  return lineageSnpData;
}
export function loadCladeSnp() {
  return cladeSnpData;
}

export function getReferenceSequence() {
  return refSeq['ref_seq'];
}

let coolColorInd = 0;
let warmColorInd = 0;
const _getLineageColor = _.memoize((group) => {
  let color;
  if (group.charAt(0) === 'A') {
    color = warmColors[warmColorInd++];

    if (warmColorInd === warmColors.length) {
      warmColorInd = 0;
    }
  } else if (group.charAt(0) === 'B') {
    color = coolColors[coolColorInd++];

    if (coolColorInd === coolColors.length) {
      coolColorInd = 0;
    }
  }
  return color;
});

const lineageColorMap = {};
lineageColorMap[GROUPS.OTHER_GROUP] = '#AAA';
lineageSnpData.forEach((lineageObj) => {
  lineageColorMap[lineageObj.lineage] = _getLineageColor(lineageObj.lineage);
});

export const getLineageColor = (lineage) => {
  return lineageColorMap[lineage];
};

let cladeColorInd = 0;
const _getCladeColor = _.memoize(() => {
  const color = cladeColorArray[cladeColorInd++];
  // If we're at the end, then loop back to the beginning
  if (cladeColorInd === cladeColorArray.length) {
    cladeColorInd = 0;
  }

  return color;
});

const cladeColorMap = {};
cladeColorMap[GROUPS.OTHER_GROUP] = '#AAA';
cladeSnpData.forEach((cladeObj) => {
  cladeColorMap[cladeObj.clade] = _getCladeColor(cladeObj.clade);
});

export const getCladeColor = (clade) => {
  return cladeColorMap[clade];
};

function getGroup(groupKey, group) {
  let snpData;
  let snpDataKey;
  if (groupKey === GROUP_KEYS.GROUP_LINEAGE) {
    snpData = lineageSnpData;
    snpDataKey = 'lineage';
  } else if (groupKey === GROUP_KEYS.GROUP_CLADE) {
    snpData = cladeSnpData;
    snpDataKey = 'clade';
  }
  let findObj = {};
  findObj[snpDataKey] = group;

  return _.findWhere(snpData, findObj);
}

export function getDnaSnpsFromGroup(groupKey, group) {
  let groupObj = getGroup(groupKey, group);
  if (groupObj === undefined) {
    return [];
  }

  let snpIds = groupObj.dna_snp_ids;
  snpIds = _.reject(snpIds, (snpId) => snpId === '');
  return _.map(snpIds, (snpId) => intToDnaSnv(parseInt(snpId)));
}

export function getGeneAaSnpsFromGroup(groupKey, group) {
  let groupObj = getGroup(groupKey, group);
  if (groupObj === undefined) {
    return [];
  }
  let snpIds = groupObj.gene_aa_snp_ids;
  snpIds = _.reject(snpIds, (snpId) => snpId === '');
  return _.map(snpIds, (snpId) => intToGeneAaSnv(parseInt(snpId)));
}

export function getProteinAaSnpsFromGroup(groupKey, group) {
  let groupObj = getGroup(groupKey, group);
  if (groupObj === undefined) {
    return [];
  }
  let snpIds = groupObj.protein_aa_snp_ids;
  snpIds = _.reject(snpIds, (snpId) => snpId === '');
  return _.map(snpIds, (snpId) => intToProteinAaSnv(parseInt(snpId)));
}

import lineageSnpData from '../../data/lineage_snp.json';
import cladeSnpData from '../../data/clade_snp.json';
import refSeq from '../../static_data/reference.json';

import _ from 'underscore';
import { intToDnaSnp, intToGeneAaSnp, intToProteinAaSnp } from './snpData';

import { GROUP_KEYS } from '../constants/config';

export function loadLineageSnp() {
  return lineageSnpData;
}
export function loadCladeSnp() {
  return cladeSnpData;
}

export function getReferenceSequence() {
  return refSeq['ref_seq'];
}

function getGroup(groupKey, group) {
  let snpData;
  if (groupKey === GROUP_KEYS.GROUP_LINEAGE) {
    snpData = lineageSnpData;
  } else if (groupKey === GROUP_KEYS.GROUP_CLADE) {
    snpData = cladeSnpData;
  }
  let findObj = {};
  findObj[groupKey] = group;

  return _.findWhere(snpData, findObj);
}

export function getDnaSnpsFromGroup(groupKey, group) {
  let groupObj = getGroup(groupKey, group);
  if (groupObj === undefined) {
    return [];
  }

  let snpIds = groupObj.dna_snp_ids;
  snpIds = _.reject(snpIds, (snpId) => snpId === '');
  return _.map(snpIds, (snpId) => intToDnaSnp(parseInt(snpId)));
}

export function getGeneAaSnpsFromGroup(groupKey, group) {
  let groupObj = getGroup(groupKey, group);
  if (groupObj === undefined) {
    return [];
  }
  let snpIds = groupObj.gene_aa_snp_ids;
  snpIds = _.reject(snpIds, (snpId) => snpId === '');
  return _.map(snpIds, (snpId) => intToGeneAaSnp(parseInt(snpId)));
}

export function getProteinAaSnpsFromGroup(groupKey, group) {
  let groupObj = getGroup(groupKey, group);
  if (groupObj === undefined) {
    return [];
  }
  let snpIds = groupObj.protein_aa_snp_ids;
  snpIds = _.reject(snpIds, (snpId) => snpId === '');
  return _.map(snpIds, (snpId) => intToProteinAaSnp(parseInt(snpId)));
}

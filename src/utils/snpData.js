/*
 * Load SNP data, and map integers -> SNP strings
 */

// import dnaSnvMap from '../../data/dna_snp_map.json';
// import geneAaSnvMap from '../../data/gene_aa_snp_map.json';
// import proteinAaSnvMap from '../../data/protein_aa_snp_map.json';

const dnaSnvMap = {};
const geneAaSnvMap = {};
const proteinAaSnvMap = {};

import _ from 'underscore';

import { getGene, getProtein } from './gene_protein';

import { snpColorArray } from '../constants/colors';
import { GROUPS } from '../constants/groups';
import { DNA_OR_AA } from '../constants/config';

class SnpDataStore {
  intToDnaSnvMap = {};
  intToGeneAaSnvMap = {};
  intToProteinAaSnvMap = {};
  snvColorMap = {};

  init() {
    //snv -> color map
    let snvColorInd = 0;
    const _getSnvColor = _.memoize(() => {
      const color = snpColorArray[snvColorInd++];
      if (snvColorInd === snpColorArray.length) {
        snvColorInd = 0;
      }
      return color;
    });

    this.snvColorMap[GROUPS.REFERENCE_GROUP] = _getSnvColor('Reference');
    this.snvColorMap[GROUPS.OTHER_GROUP] = '#AAA';
    this.snvColorMap[GROUPS.NONE_GROUP] = '#AAA';
    this.snvColorMap[GROUPS.ALL_OTHER_GROUP] = '#AAA';
    Object.keys(dnaSnvMap).forEach((snv) => {
      this.snvColorMap[snv] = _getSnvColor(snv);
    });
    Object.keys(geneAaSnvMap).forEach((snv) => {
      this.snvColorMap[snv] = _getSnvColor(snv);
    });
    Object.keys(proteinAaSnvMap).forEach((snv) => {
      this.snvColorMap[snv] = _getSnvColor(snv);
    });

    //remap so its integer -> SNP
    let snvId, split, aaRangeInd;

    Object.keys(dnaSnvMap).forEach((snv) => {
      snvId = parseInt(dnaSnvMap[snv]);
      this.intToDnaSnvMap[snvId] = {};
      // Store the entire SNV string
      this.intToDnaSnvMap[snvId]['snp_str'] = snv;

      // Each SNV is broken up by pos|ref|alt
      split = snv.split('|');

      // Store all the parts
      // Positions are 1-indexed
      this.intToDnaSnvMap[snvId]['pos'] = parseInt(split[0]);
      this.intToDnaSnvMap[snvId]['ref'] = split[1];
      this.intToDnaSnvMap[snvId]['alt'] = split[2];
    });
    Object.keys(geneAaSnvMap).forEach((snv) => {
      snvId = parseInt(geneAaSnvMap[snv]);
      this.intToGeneAaSnvMap[snvId] = {};
      // Store the entire SNV string
      this.intToGeneAaSnvMap[snvId]['snp_str'] = snv;

      // Each SNV is broken up by gene|pos|ref|alt
      split = snv.split('|');

      // Store all the parts
      this.intToGeneAaSnvMap[snvId]['gene'] = split[0];
      this.intToGeneAaSnvMap[snvId]['pos'] = parseInt(split[1]);
      this.intToGeneAaSnvMap[snvId]['ref'] = split[2];
      this.intToGeneAaSnvMap[snvId]['alt'] = split[3];

      // Get coordinates in NT (from start of codon)
      aaRangeInd = getGene(split[0]).aa_ranges.reduce(
        (_aaRangeInd, range, ind) => {
          return this.intToGeneAaSnvMap[snvId]['pos'] >= range[0] &&
            this.intToGeneAaSnvMap[snvId]['pos'] <= range[1]
            ? ind
            : _aaRangeInd;
        },
        0
      );
      this.intToGeneAaSnvMap[snvId]['nt_pos'] =
        getGene(split[0]).ranges[aaRangeInd][0] +
        (this.intToGeneAaSnvMap[snvId]['pos'] -
          getGene(split[0]).aa_ranges[aaRangeInd][0]) *
          3;
    });
    Object.keys(proteinAaSnvMap).forEach((snv) => {
      snvId = parseInt(proteinAaSnvMap[snv]);
      this.intToProteinAaSnvMap[snvId] = {};
      // Store the entire SNV string
      this.intToProteinAaSnvMap[snvId]['snp_str'] = snv;

      // Each SNV is broken up by gene|pos|ref|alt
      split = snv.split('|');

      // Store all the parts
      this.intToProteinAaSnvMap[snvId]['protein'] = split[0];
      this.intToProteinAaSnvMap[snvId]['pos'] = parseInt(split[1]);
      this.intToProteinAaSnvMap[snvId]['ref'] = split[2];
      this.intToProteinAaSnvMap[snvId]['alt'] = split[3];

      // Get coordinates in NT (from start of codon)
      aaRangeInd = getProtein(split[0]).aa_ranges.reduce(
        (_aaRangeInd, range, ind) => {
          return this.intToProteinAaSnvMap[snvId]['pos'] >= range[0] &&
            this.intToProteinAaSnvMap[snvId]['pos'] <= range[1]
            ? ind
            : _aaRangeInd;
        },
        0
      );
      this.intToProteinAaSnvMap[snvId]['nt_pos'] =
        getProtein(split[0]).ranges[aaRangeInd][0] +
        (this.intToProteinAaSnvMap[snvId]['pos'] -
          getProtein(split[0]).aa_ranges[aaRangeInd][0]) *
          3;
    });
  }

  getSnvColor(snv) {
    return this.snvColorMap[snv];
  }
  intToDnaSnv(dnaSnvId) {
    if (dnaSnvId === -1) {
      return { snp_str: GROUPS.REFERENCE_GROUP };
    }
    return this.intToDnaSnvMap[dnaSnvId];
  }
  intToGeneAaSnv(aaSnvId) {
    if (aaSnvId === -1) {
      return { snp_str: GROUPS.REFERENCE_GROUP };
    }
    return this.intToGeneAaSnvMap[aaSnvId];
  }
  intToProteinAaSnv(aaSnvId) {
    if (aaSnvId === -1) {
      return { snp_str: GROUPS.REFERENCE_GROUP };
    }
    return this.intToProteinAaSnvMap[aaSnvId];
  }

  dnaSnvToInt(dnaSnv) {
    return dnaSnvMap[dnaSnv];
  }
  geneAaSnvToInt(geneAaSnv) {
    return geneAaSnvMap[geneAaSnv];
  }
  proteinAaSnvToInt(proteinAaSnv) {
    return proteinAaSnvMap[proteinAaSnv];
  }

  formatSnv(snvStr, dnaOrAa) {
    // Don't do this if it's a special group
    if (Object.values(GROUPS).includes(snvStr)) {
      return snvStr;
    }

    // Print as REF POS ALT
    // i.e., 23403|A|G -> A23403G, S|614|D|G -> S · D614G
    const chunks = snvStr.split('|');
    if (dnaOrAa === DNA_OR_AA.DNA) {
      return `${chunks[1]}${chunks[0]}${chunks[2]}`;
    } else if (dnaOrAa === DNA_OR_AA.AA) {
      return `${chunks[0]} · ${chunks[2]}${chunks[1]}${chunks[3]}`;
    }
  }
}

export const snpDataStoreInstance = new SnpDataStore();

/*
 * Load SNP data, and map integers -> SNP strings
 */

// import dnaSnvMap from '../../data/dna_snp_map.json';
// import geneAaSnvMap from '../../data/gene_aa_snp_map.json';
// import proteinAaSnvMap from '../../data/protein_aa_snp_map.json';

import _ from 'underscore';

import { getGene, getProtein } from '../utils/gene_protein';

import { snpColorArray } from '../constants/colors';
import { GROUPS } from '../constants/groups';
import { asyncDataStoreInstance } from '../components/App';

export class SnpDataStore {
  intToDnaSnvMap = {
    '-1': { snp_str: GROUPS.REFERENCE_GROUP },
  };
  intToGeneAaSnvMap = {
    '-1': { snp_str: GROUPS.REFERENCE_GROUP },
  };
  intToProteinAaSnvMap = {
    '-1': { snp_str: GROUPS.REFERENCE_GROUP },
  };
  snvColorMap = {};
  dnaSnvMap = {};
  geneAaSnvMap = {};
  proteinAaSnvMap = {};

  constructor() {
    this.dnaSnvMap = asyncDataStoreInstance.data.dna_snp_map;
    this.geneAaSnvMap = asyncDataStoreInstance.data.gene_aa_snp_map;
    this.proteinAaSnvMap = asyncDataStoreInstance.data.protein_aa_snp_map;
    // debugger;

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
    Object.keys(this.dnaSnvMap).forEach((snv) => {
      this.snvColorMap[snv] = _getSnvColor(snv);
    });
    Object.keys(this.geneAaSnvMap).forEach((snv) => {
      this.snvColorMap[snv] = _getSnvColor(snv);
    });
    Object.keys(this.proteinAaSnvMap).forEach((snv) => {
      this.snvColorMap[snv] = _getSnvColor(snv);
    });

    //remap so its integer -> SNP
    let snvId, split, aaRangeInd;

    Object.keys(this.dnaSnvMap).forEach((snv) => {
      snvId = parseInt(this.dnaSnvMap[snv]);
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
    Object.keys(this.geneAaSnvMap).forEach((snv) => {
      snvId = parseInt(this.geneAaSnvMap[snv]);
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
    Object.keys(this.proteinAaSnvMap).forEach((snv) => {
      snvId = parseInt(this.proteinAaSnvMap[snv]);
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
    return this.intToDnaSnvMap[dnaSnvId];
  }
  intToGeneAaSnv(aaSnvId) {
    return this.intToGeneAaSnvMap[aaSnvId];
  }
  intToProteinAaSnv(aaSnvId) {
    return this.intToProteinAaSnvMap[aaSnvId];
  }

  dnaSnvToInt(dnaSnv) {
    return this.dnaSnvMap[dnaSnv];
  }
  geneAaSnvToInt(geneAaSnv) {
    return this.geneAaSnvMap[geneAaSnv];
  }
  proteinAaSnvToInt(proteinAaSnv) {
    return this.proteinAaSnvMap[proteinAaSnv];
  }
}

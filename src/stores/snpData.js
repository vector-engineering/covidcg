/*
 * Load SNP data, and map integers -> SNP strings
 */

import _ from 'underscore';

import { getGene, getProtein } from '../utils/gene_protein';
import { formatSnv } from '../utils/snpUtils';
import { snpColorArray } from '../constants/colors';
import { GROUPS, DNA_OR_AA, COORDINATE_MODES } from '../constants/defs.json';
import { asyncDataStoreInstance } from '../components/App';

export class SnpDataStore {
  intToDnaSnvMap;
  intToGeneAaSnvMap;
  intToProteinAaSnvMap;
  dnaSnvMap;
  geneAaSnvMap;
  proteinAaSnvMap;
  snvColorMap;

  constructor() {
    this.intToDnaSnvMap = {
      '-1': { snp_str: GROUPS.REFERENCE_GROUP },
    };
    this.intToGeneAaSnvMap = {
      '-1': { snp_str: GROUPS.REFERENCE_GROUP },
    };
    this.intToProteinAaSnvMap = {
      '-1': { snp_str: GROUPS.REFERENCE_GROUP },
    };
    this.dnaSnvMap = {};
    this.geneAaSnvMap = {};
    this.proteinAaSnvMap = {};

    this.snvColorMap = {};
  }

  init() {
    this.dnaSnvMap = asyncDataStoreInstance.data.metadata_map.dna_snp;
    this.geneAaSnvMap = asyncDataStoreInstance.data.metadata_map.gene_aa_snp;
    this.proteinAaSnvMap =
      asyncDataStoreInstance.data.metadata_map.protein_aa_snp;
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
      this.intToDnaSnvMap[snvId]['name'] = formatSnv(snv, DNA_OR_AA.DNA);
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
      this.intToGeneAaSnvMap[snvId]['name'] = formatSnv(snv, DNA_OR_AA.AA);

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
        getGene(split[0]).segments[aaRangeInd][0] +
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
      this.intToProteinAaSnvMap[snvId]['name'] = formatSnv(snv, DNA_OR_AA.AA);

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
        getProtein(split[0]).segments[aaRangeInd][0] +
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
  intToSnv(dnaOrAa, coordinateMode, snvId) {
    if (dnaOrAa === DNA_OR_AA.DNA) {
      return this.intToDnaSnv(snvId);
    } else {
      if (coordinateMode === COORDINATE_MODES.COORD_GENE) {
        return this.intToGeneAaSnv(snvId);
      } else {
        return this.intToProteinAaSnv(snvId);
      }
    }
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
  snvToInt(dnaOrAa, coordinateMode, snv) {
    if (dnaOrAa === DNA_OR_AA.DNA) {
      return this.dnaSnvToInt(snv);
    } else {
      if (coordinateMode === COORDINATE_MODES.COORD_GENE) {
        return this.geneAaSnvToInt(snv);
      } else {
        return this.proteinAaSnvToInt(snv);
      }
    }
  }
}

/*
 * Load mutation data, and map integers -> mutation strings
 */

import { config } from '../config';
import { getGene, getProtein } from '../utils/gene_protein';
import { formatMutation } from '../utils/mutationUtils';
import { memoize } from '../utils/func';
import { mutationColorArray } from '../constants/colors';
import { GROUPS, DNA_OR_AA, COORDINATE_MODES } from '../constants/defs.json';
import { asyncDataStoreInstance } from '../components/App';
import { rootStoreInstance } from './rootStore';

export class MutationDataStore {
  intToDnaMutationMap;
  intToGeneAaMutationMap;
  intToProteinAaMutationMap;
  dnaMutationMap;
  geneAaMutationMap;
  proteinAaMutationMap;
  mutationColorMap;

  constructor() {
    this.intToDnaMutationMap = {
      '-1': { mutation_str: GROUPS.REFERENCE_GROUP },
    };
    this.intToGeneAaMutationMap = {
      '-1': { mutation_str: GROUPS.REFERENCE_GROUP },
    };
    this.intToProteinAaMutationMap = {
      '-1': { mutation_str: GROUPS.REFERENCE_GROUP },
    };
    this.dnaMutationMap = {};
    this.geneAaMutationMap = {};
    this.proteinAaMutationMap = {};

    this.mutationColorMap = {};
  }

  init() {
    this.dnaMutationMap = asyncDataStoreInstance.data.metadata_map.dna_mutation;
    this.geneAaMutationMap =
      asyncDataStoreInstance.data.metadata_map.gene_aa_mutation;
    this.proteinAaMutationMap =
      asyncDataStoreInstance.data.metadata_map.protein_aa_mutation;
    // debugger;

    //mutation -> color map
    let mutationColorInd = 0;
    const _getMutationColor = memoize(() => {
      const color = mutationColorArray[mutationColorInd++];
      if (mutationColorInd === mutationColorArray.length) {
        mutationColorInd = 0;
      }
      return color;
    });

    this.mutationColorMap[GROUPS.REFERENCE_GROUP] =
      _getMutationColor('Reference');
    this.mutationColorMap[GROUPS.OTHER_GROUP] = '#AAA';
    this.mutationColorMap[GROUPS.NONE_GROUP] = '#AAA';
    this.mutationColorMap[GROUPS.ALL_OTHER_GROUP] = '#AAA';
    Object.keys(this.dnaMutationMap).forEach((mut) => {
      this.mutationColorMap[mut] = _getMutationColor(mut);
    });
    Object.keys(this.geneAaMutationMap).forEach((mut) => {
      this.mutationColorMap[mut] = _getMutationColor(mut);
    });
    Object.keys(this.proteinAaMutationMap).forEach((mut) => {
      this.mutationColorMap[mut] = _getMutationColor(mut);
    });

    //remap so its integer -> mutation
    let mutationId, split, aaRangeInd;

    Object.keys(this.dnaMutationMap).forEach((mut) => {
      mutationId = parseInt(this.dnaMutationMap[mut]);
      this.intToDnaMutationMap[mutationId] = {};
      // Store the entire mutation string
      this.intToDnaMutationMap[mutationId]['mutation_str'] = mut;

      // Each mut is broken up by pos|ref|alt
      split = mut.split('|');

      // Store all the parts
      // Positions are 1-indexed
      this.intToDnaMutationMap[mutationId]['pos'] = parseInt(split[0]);
      this.intToDnaMutationMap[mutationId]['ref'] = split[1];
      this.intToDnaMutationMap[mutationId]['alt'] = split[2];
      this.intToDnaMutationMap[mutationId]['name'] = formatMutation(
        mut,
        DNA_OR_AA.DNA
      );
      this.intToDnaMutationMap[mutationId]['color'] =
        this.mutationColorMap[mut];
    });
    Object.keys(this.geneAaMutationMap).forEach((mut) => {
      mutationId = parseInt(this.geneAaMutationMap[mut]);
      this.intToGeneAaMutationMap[mutationId] = {};
      // Store the entire mutation string
      this.intToGeneAaMutationMap[mutationId]['mutation_str'] = mut;

      // Each mutation is broken up by gene|pos|ref|alt
      split = mut.split('|');

      // Store all the parts
      this.intToGeneAaMutationMap[mutationId]['gene'] = split[0];
      this.intToGeneAaMutationMap[mutationId]['pos'] = parseInt(split[1]);
      this.intToGeneAaMutationMap[mutationId]['ref'] = split[2];
      this.intToGeneAaMutationMap[mutationId]['alt'] = split[3];
      this.intToGeneAaMutationMap[mutationId]['name'] = formatMutation(
        mut,
        DNA_OR_AA.AA
      );
      this.intToGeneAaMutationMap[mutationId]['color'] =
        this.mutationColorMap[mut];

      // Get coordinates in NT (from start of codon)
      let gene;
      if (config.virus === 'sars2') {
        gene = getGene(split[0]);
      } else if (config.virus === 'rsv') {
        gene = getGene(
          split[0],
          rootStoreInstance.configStore.selectedReference
        );
      }

      aaRangeInd = gene.aa_ranges.reduce((_aaRangeInd, range, ind) => {
        return this.intToGeneAaMutationMap[mutationId]['pos'] >= range[0] &&
          this.intToGeneAaMutationMap[mutationId]['pos'] <= range[1]
          ? ind
          : _aaRangeInd;
      }, 0);

      this.intToGeneAaMutationMap[mutationId]['nt_pos'] =
        gene.segments[aaRangeInd][0] +
        (this.intToGeneAaMutationMap[mutationId]['pos'] -
          gene.aa_ranges[aaRangeInd][0]) *
          3;
    });
    Object.keys(this.proteinAaMutationMap).forEach((mut) => {
      mutationId = parseInt(this.proteinAaMutationMap[mut]);
      this.intToProteinAaMutationMap[mutationId] = {};
      // Store the entire mutation string
      this.intToProteinAaMutationMap[mutationId]['mutation_str'] = mut;

      // Each mutation is broken up by gene|pos|ref|alt
      split = mut.split('|');

      // Store all the parts
      this.intToProteinAaMutationMap[mutationId]['protein'] = split[0];
      this.intToProteinAaMutationMap[mutationId]['pos'] = parseInt(split[1]);
      this.intToProteinAaMutationMap[mutationId]['ref'] = split[2];
      this.intToProteinAaMutationMap[mutationId]['alt'] = split[3];
      this.intToProteinAaMutationMap[mutationId]['name'] = formatMutation(
        mut,
        DNA_OR_AA.AA
      );
      this.intToProteinAaMutationMap[mutationId]['color'] =
        this.mutationColorMap[mut];

      // Get coordinates in NT (from start of codon)

      let protein;
      if (config.virus === 'sars2') {
        protein = getProtein(split[0]);
      } else if (config.virus === 'rsv') {
        protein = getProtein(
          split[0],
          rootStoreInstance.configStore.selectedReference
        );
      }

      aaRangeInd = protein.aa_ranges.reduce((_aaRangeInd, range, ind) => {
        return this.intToProteinAaMutationMap[mutationId]['pos'] >= range[0] &&
          this.intToProteinAaMutationMap[mutationId]['pos'] <= range[1]
          ? ind
          : _aaRangeInd;
      }, 0);
      this.intToProteinAaMutationMap[mutationId]['nt_pos'] =
        protein.segments[aaRangeInd][0] +
        (this.intToProteinAaMutationMap[mutationId]['pos'] -
          protein.aa_ranges[aaRangeInd][0]) *
          3;
    });
  }

  getMutationColor(mut) {
    return this.mutationColorMap[mut];
  }
  intToDnaMutation(dnaMutationId) {
    return this.intToDnaMutationMap[dnaMutationId];
  }
  intToGeneAaMutation(aaMutationId) {
    return this.intToGeneAaMutationMap[aaMutationId];
  }
  intToProteinAaMutation(aaMutationId) {
    return this.intToProteinAaMutationMap[aaMutationId];
  }
  intToMutation(dnaOrAa, coordinateMode, mutationId) {
    if (dnaOrAa === DNA_OR_AA.DNA) {
      return this.intToDnaMutation(mutationId);
    } else {
      if (coordinateMode === COORDINATE_MODES.COORD_GENE) {
        return this.intToGeneAaMutation(mutationId);
      } else {
        return this.intToProteinAaMutation(mutationId);
      }
    }
  }

  dnaMutationToInt(dnaMutation) {
    return this.dnaMutationMap[dnaMutation];
  }
  geneAaMutationToInt(geneAaMutation) {
    return this.geneAaMutationMap[geneAaMutation];
  }
  proteinAaMutationToInt(proteinAaMutation) {
    return this.proteinAaMutationMap[proteinAaMutation];
  }
  mutationToInt(dnaOrAa, coordinateMode, mut) {
    if (dnaOrAa === DNA_OR_AA.DNA) {
      return this.dnaMutationToInt(mut);
    } else {
      if (coordinateMode === COORDINATE_MODES.COORD_GENE) {
        return this.geneAaMutationToInt(mut);
      } else {
        return this.proteinAaMutationToInt(mut);
      }
    }
  }
}

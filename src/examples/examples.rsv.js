import { getLocationByNameAndLevel } from '../utils/location';
import { getGene } from '../utils/gene_protein';
import { todayISO } from '../utils/date';

import {
  GROUP_MUTATION,
  DNA_OR_AA,
  COORDINATE_MODES,
  NORM_MODES,
  COUNT_MODES,
  DATE_BINS,
  TABS,
} from '../constants/defs.json';

import { config } from '../config';

export default function examples({ selectTree }) {
  return [
    //   {
    //     title: 'Global Subtypes',
    //     description: 'View the growth of subtypes over all locations',
    //     // image:
    //     //   'https://storage.googleapis.com/ve-public/example/global_lineages.png',
    //     settings: {
    //       plotSettings: {
    //         groupStackNormMode: NORM_MODES.NORM_PERCENTAGES,
    //         groupStackCountMode: COUNT_MODES.COUNT_NEW,
    //         groupStackDateBin: DATE_BINS.DATE_BIN_YEAR,
    //       },
    //       UI: {
    //         activeTab: TABS.TAB_COMPARE_GROUPS,
    //       },
    //       config: {
    //         groupKey: 'subtype',
    //         dnaOrAa: DNA_OR_AA.DNA,
    //         coordinateMode: COORDINATE_MODES.COORD_GENE,
    //         selectedLocationNodes: selectTree.children, // select root
    //         startDate: config.min_date,
    //         endDate: todayISO(),
    //         selectedGroupFields: {},
    //       },
    //     },
    //   },
    {
      title: 'Global Genotypes',
      description: 'View the growth of genotypes over all locations',
      // image:
      //   'https://storage.googleapis.com/ve-public/example/global_lineages.png',
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_PERCENTAGES,
          groupStackCountMode: COUNT_MODES.COUNT_NEW,
          groupStackDateBin: DATE_BINS.DATE_BIN_YEAR,
        },
        UI: {
          activeTab: TABS.TAB_COMPARE_GROUPS,
        },
        config: {
          groupKey: config.group_cols.genotype.name,
          dnaOrAa: DNA_OR_AA.DNA,
          coordinateMode: COORDINATE_MODES.COORD_GENE,
          selectedLocationNodes: selectTree.children, // select root
          startDate: config.min_date,
          endDate: todayISO(),
          selectedGroupFields: {},
        },
      },
    },
    {
      title: 'Europe: Fusion A',
      description:
        'View mutations in the fusion (F) gene in A subtypes, in Europe over all time',
      // image:
      //   'https://storage.googleapis.com/ve-public/example/global_lineages.png',
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_PERCENTAGES,
          groupStackCountMode: COUNT_MODES.COUNT_NEW,
          groupStackDateBin: DATE_BINS.DATE_BIN_YEAR,
        },
        UI: {
          activeTab: TABS.TAB_COMPARE_GROUPS,
        },
        config: {
          groupKey: GROUP_MUTATION,
          dnaOrAa: DNA_OR_AA.AA,
          coordinateMode: COORDINATE_MODES.COORD_GENE,
          selectedReference: 'NC_038235.1',
          selectedGene: getGene('F', 'NC_038235.1'),
          selectedLocationNodes: [
            getLocationByNameAndLevel(selectTree, 'Europe', 'region')[0],
          ],
          startDate: config.min_date,
          endDate: todayISO(),
          selectedGroupFields: {
            subtype: ['A'],
          },
        },
      },
    },
    {
      title: 'North America: Fusion B',
      description:
        'View mutations in the fusion (F) gene in B subtypes, in North America over all time',
      // image:
      //   'https://storage.googleapis.com/ve-public/example/global_lineages.png',
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_PERCENTAGES,
          groupStackCountMode: COUNT_MODES.COUNT_NEW,
          groupStackDateBin: DATE_BINS.DATE_BIN_YEAR,
        },
        UI: {
          activeTab: TABS.TAB_COMPARE_GROUPS,
        },
        config: {
          groupKey: GROUP_MUTATION,
          dnaOrAa: DNA_OR_AA.AA,
          coordinateMode: COORDINATE_MODES.COORD_GENE,
          selectedReference: 'NC_001781.1',
          selectedGene: getGene('F', 'NC_001781.1'),
          selectedLocationNodes: [
            getLocationByNameAndLevel(selectTree, 'North America', 'region')[0],
          ],
          startDate: config.min_date,
          endDate: todayISO(),
          selectedGroupFields: {
            subtype: ['B'],
          },
        },
      },
    },
  ];
}

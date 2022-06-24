import { getGene } from '../utils/gene_protein';
import { getLocationByNameAndLevel } from '../utils/location';
import { queryPrimers } from '../utils/primer';
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
    {
      title: 'Global Lineages',
      description: 'View the growth of the B lineage family over all locations',
      image:
        'https://storage.googleapis.com/ve-public/example/global_lineages.png',
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_PERCENTAGES,
          groupStackCountMode: COUNT_MODES.COUNT_NEW,
          groupStackDateBin: DATE_BINS.DATE_BIN_WEEK,
        },
        UI: {
          activeTab: TABS.TAB_COMPARE_GROUPS,
        },
        config: {
          groupKey: config.group_cols.lineage.name,
          dnaOrAa: DNA_OR_AA.DNA,
          selectedGene: getGene('S'),
          coordinateMode: COORDINATE_MODES.COORD_GENE,
          selectedLocationNodes: [selectTree], // select root
          startDate: config.min_date,
          endDate: todayISO(),
        },
      },
    },
    // {
    //   title: 'Lineages in China - Beijing Xinfadi Market',
    //   description:
    //     "New lineages in uncovered in early June in Beijing's Xinfadi Market may have been circulating in China in March",
    //   image: 'https://storage.googleapis.com/ve-public/example/xinfadi.png',
    //   settings: {
    //     plotSettings: {
    //       groupStackNormMode: NORM_MODES.NORM_COUNTS,
    //       groupStackCountMode: COUNT_MODES.COUNT_NEW,
    //       groupStackDateBin: DATE_BINS.DATE_BIN_DAY,
    //     },
    //     UI: {
    //       activeTab: TABS.TAB_COMPARE_GROUPS,
    //     },
    //     config: {
    //       groupKey: config.group_cols.lineage.name,
    //       dnaOrAa: DNA_OR_AA.AA,
    //       selectedGene: getGene('S'),
    //       coordinateMode: COORDINATE_MODES.COORD_GENE,
    //       selectedLocationNodes: [getLocationByNameAndLevel(selectTree, 'China', 'country')[0]],
    //     },
    //   },
    // },
    {
      title: 'Lineages in Iceland',
      description:
        'Lineages sequenced in Iceland during the early stages of the pandemic',
      image: 'https://storage.googleapis.com/ve-public/example/iceland.png',
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_COUNTS,
          groupStackCountMode: COUNT_MODES.COUNT_NEW,
          groupStackDateBin: DATE_BINS.DATE_BIN_DAY,
        },
        UI: {
          activeTab: TABS.TAB_COMPARE_GROUPS,
        },
        config: {
          groupKey: config.group_cols.lineage.name,
          dnaOrAa: DNA_OR_AA.AA,
          selectedGene: getGene('S'),
          coordinateMode: COORDINATE_MODES.COORD_GENE,
          selectedLocationNodes: [
            getLocationByNameAndLevel(selectTree, 'Iceland', 'country')[0],
          ],
          startDate: config.min_date,
          endDate: todayISO(),
        },
      },
    },
    {
      title: 'Rise of Spike D614G mutation in West Coast USA',
      description:
        'The Spike D614G mutation has accumulated in frequency in the West Coast states of the USA',
      image:
        'https://storage.googleapis.com/ve-public/example/d614g_west_coast.png',
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_PERCENTAGES,
          groupStackCountMode: COUNT_MODES.COUNT_NEW,
          groupStackDateBin: DATE_BINS.DATE_BIN_DAY,
        },
        UI: {
          activeTab: TABS.TAB_COMPARE_GROUPS,
        },
        config: {
          groupKey: GROUP_MUTATION,
          dnaOrAa: DNA_OR_AA.AA,
          selectedGene: getGene('S'),
          coordinateMode: COORDINATE_MODES.COORD_GENE,
          selectedLocationNodes: [
            getLocationByNameAndLevel(selectTree, 'Washington', 'division')[0],
            getLocationByNameAndLevel(selectTree, 'Oregon', 'division')[0],
            getLocationByNameAndLevel(selectTree, 'California', 'division')[0],
          ],
          startDate: '2020-03-01',
          endDate: '2020-06-01',
          selectedGroups: [{ group: 'S|614|D|G' }],
        },
      },
    },
    {
      title: 'Prevalence of Spike D614G in various US States',
      description:
        'The proportion of sequences with the Spike D614G mutation varies between US States',
      image:
        'https://storage.googleapis.com/ve-public/example/d614g_us_states.png',
      settings: {
        plotSettings: {
          locationDateNormMode: NORM_MODES.NORM_PERCENTAGES,
          locationDateCountMode: COUNT_MODES.COUNT_CUMULATIVE,
          locationDateDateBin: DATE_BINS.DATE_BIN_MONTH,
        },
        UI: {
          activeTab: TABS.TAB_COMPARE_LOCATIONS,
        },
        config: {
          groupKey: GROUP_MUTATION,
          dnaOrAa: DNA_OR_AA.AA,
          selectedGene: getGene('S'),
          coordinateMode: COORDINATE_MODES.COORD_GENE,
          selectedGroups: [{ group: 'S|614|D|G' }],
          selectedLocationNodes: [
            getLocationByNameAndLevel(selectTree, 'Washington', 'division')[0],
            getLocationByNameAndLevel(selectTree, 'California', 'division')[0],
            getLocationByNameAndLevel(selectTree, 'Florida', 'division')[0],
            getLocationByNameAndLevel(
              selectTree,
              'Massachusetts',
              'division'
            )[0],
            getLocationByNameAndLevel(selectTree, 'Michigan', 'division')[0],
            getLocationByNameAndLevel(selectTree, 'New York', 'division')[0],
            getLocationByNameAndLevel(selectTree, 'Texas', 'division')[0],
            getLocationByNameAndLevel(selectTree, 'Wisconsin', 'division')[0],
          ],
          startDate: config.min_date,
          endDate: todayISO(),
        },
      },
    },
    {
      title: 'Diagnostics: NT mutations in US CDC qPCR primer/probe sequences',
      description:
        'Prevalence of any mutations present within the US CDC primer and probe sequences (N1 + N2), for sequences in the US',
      image:
        'https://storage.googleapis.com/ve-public/example/us_cdc_primer.png',
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_PERCENTAGES,
          groupStackCountMode: COUNT_MODES.COUNT_CUMULATIVE,
          groupStackDateBin: DATE_BINS.DATE_BIN_WEEK,
        },
        UI: {
          activeTab: TABS.TAB_COMPARE_GROUPS,
        },
        config: {
          groupKey: GROUP_MUTATION,
          dnaOrAa: DNA_OR_AA.DNA,
          coordinateMode: COORDINATE_MODES.COORD_PRIMER,
          selectedLocationNodes: [
            getLocationByNameAndLevel(selectTree, 'USA', 'country')[0],
          ],
          selectedPrimers: []
            .concat(
              queryPrimers({
                Institution: 'US CDC',
                Name: '2019-nCoV-N1-F',
              })
            )
            .concat(
              queryPrimers({
                Institution: 'US CDC',
                Name: '2019-nCoV-N1-R',
              })
            )
            .concat(
              queryPrimers({
                Institution: 'US CDC',
                Name: '2019-nCoV-N1-P',
              })
            )
            .concat(
              queryPrimers({
                Institution: 'US CDC',
                Name: '2019-nCoV-N2-F',
              })
            )
            .concat(
              queryPrimers({
                Institution: 'US CDC',
                Name: '2019-nCoV-N2-R',
              })
            )
            .concat(
              queryPrimers({
                Institution: 'US CDC',
                Name: '2019-nCoV-N2-P',
              })
            ),
          startDate: config.min_date,
          endDate: todayISO(),
        },
      },
    },
    {
      title: 'Co-occurrence of R203K and G204R in N gene',
      description:
        'Two mutations in the N gene, R203K and G204R, co-occur with each other. These mutations are associated with the B.1.1 lineage.',
      image:
        'https://storage.googleapis.com/ve-public/example/n_203_204_coocurrence.png',
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_PERCENTAGES,
          groupStackCountMode: COUNT_MODES.COUNT_NEW,
          groupStackDateBin: DATE_BINS.DATE_BIN_WEEK,
        },
        UI: {
          activeTab: TABS.TAB_COMPARE_GROUPS,
        },
        config: {
          groupKey: GROUP_MUTATION,
          dnaOrAa: DNA_OR_AA.AA,
          selectedGene: getGene('N'),
          coordinateMode: COORDINATE_MODES.COORD_GENE,
          selectedLocationNodes: [
            getLocationByNameAndLevel(selectTree, 'USA', 'country')[0],
          ],
          selectedGroups: [{ group: 'N|203|R|K' }, { group: 'N|204|G|R' }],
          startDate: config.min_date,
          endDate: todayISO(),
        },
      },
    },
    {
      title: 'Emergence of Spike D614G in Europe/North America',
      description: '',
      image:
        'https://storage.googleapis.com/ve-public/example/d614g_europe_na.png',
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_PERCENTAGES,
          groupStackCountMode: COUNT_MODES.COUNT_NEW,
          groupStackDateBin: DATE_BINS.DATE_BIN_WEEK,
        },
        UI: {
          activeTab: TABS.TAB_COMPARE_GROUPS,
        },
        config: {
          groupKey: GROUP_MUTATION,
          dnaOrAa: DNA_OR_AA.AA,
          selectedGene: getGene('S'),
          coordinateMode: COORDINATE_MODES.COORD_GENE,
          selectedLocationNodes: [
            getLocationByNameAndLevel(selectTree, 'North America', 'region')[0],
            getLocationByNameAndLevel(selectTree, 'Europe', 'region')[0],
          ],
          selectedGroups: [{ group: 'S|614|D|G' }],
          startDate: config.min_date,
          endDate: '2020-12-31',
        },
      },
    },
  ];
}

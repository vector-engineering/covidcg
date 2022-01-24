import { getGene } from '../../utils/gene_protein';
import { getLocationByNameAndLevel } from '../../utils/location';
import { queryPrimers } from '../../utils/primer';
import { todayISO } from '../../utils/date';

<<<<<<< HEAD
import GlobalLineagesImage from '../../assets/analysis_screens/global_lineages.png';
// import XinfadiImage from '../../assets/analysis_screens/xinfadi.png';
import IcelandImage from '../../assets/analysis_screens/iceland.png';
import D614GWestCoastImage from '../../assets/analysis_screens/d614g_west_coast.png';
import D614GUSStatesImage from '../../assets/analysis_screens/d614g_us_states.png';
import USCDCPrimerImage from '../../assets/analysis_screens/us_cdc_primer.png';
import D614GEuropeNAImage from '../../assets/analysis_screens/d614g_europe_na.png';
import N203204Image from '../../assets/analysis_screens/n_203_204_coocurrence.png';

=======
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
import {
  GROUP_MUTATION,
  DNA_OR_AA,
  COORDINATE_MODES,
  NORM_MODES,
  COUNT_MODES,
  DATE_BINS,
  TABS,
  MIN_DATE,
} from '../../constants/defs.json';

import { config } from '../../config';

<<<<<<< HEAD
const min_date = MIN_DATE.SARS2;

export default function examples(selectTree) {
=======
export const examples = ({ selectTree }) => {
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
  return [
    {
      title: 'Global Lineages',
      description: 'View the growth of the B lineage family over all locations',
<<<<<<< HEAD
      image: GlobalLineagesImage,
=======
      image:
        'https://storage.googleapis.com/ve-public/example/global_lineages.png',
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
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
<<<<<<< HEAD
          selectedLocationNodes: [selectTree], // select root
          startDate: min_date,
=======
          selectedLocationNodes: selectTree.children, // select root
          startDate: MIN_DATE,
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
          endDate: todayISO(),
        },
      },
    },
    // {
    //   title: 'Lineages in China - Beijing Xinfadi Market',
    //   description:
    //     "New lineages in uncovered in early June in Beijing's Xinfadi Market may have been circulating in China in March",
<<<<<<< HEAD
    //   image: XinfadiImage,
=======
    //   image: 'https://storage.googleapis.com/ve-public/example/xinfadi.png',
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
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
<<<<<<< HEAD
      image: IcelandImage,
=======
      image: 'https://storage.googleapis.com/ve-public/example/iceland.png',
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
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
<<<<<<< HEAD
          startDate: min_date,
=======
          startDate: MIN_DATE,
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
          endDate: todayISO(),
        },
      },
    },
    {
      title: 'Rise of Spike D614G mutation in West Coast USA',
      description:
        'The Spike D614G mutation has accumulated in frequency in the West Coast states of the USA',
<<<<<<< HEAD
      image: D614GWestCoastImage,
=======
      image:
        'https://storage.googleapis.com/ve-public/example/d614g_west_coast.png',
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
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
<<<<<<< HEAD
      image: D614GUSStatesImage,
=======
      image:
        'https://storage.googleapis.com/ve-public/example/d614g_us_states.png',
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
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
<<<<<<< HEAD
          startDate: min_date,
=======
          startDate: MIN_DATE,
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
          endDate: todayISO(),
        },
      },
    },
    {
<<<<<<< HEAD
      title: 'Diagnostics: NT MUTATIONs in US CDC qPCR primer/probe sequences',
      description:
        'Prevalence of any mutations present within the US CDC primer and probe sequences (N1 + N2), for sequences in the US',
      image: USCDCPrimerImage,
=======
      title: 'Diagnostics: NT mutations in US CDC qPCR primer/probe sequences',
      description:
        'Prevalence of any mutations present within the US CDC primer and probe sequences (N1 + N2), for sequences in the US',
      image:
        'https://storage.googleapis.com/ve-public/example/us_cdc_primer.png',
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
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
<<<<<<< HEAD
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
          startDate: min_date,
=======
              queryPrimers({ Institution: 'US CDC', Name: '2019-nCoV-N1-F' })
            )
            .concat(
              queryPrimers({ Institution: 'US CDC', Name: '2019-nCoV-N1-R' })
            )
            .concat(
              queryPrimers({ Institution: 'US CDC', Name: '2019-nCoV-N1-P' })
            )
            .concat(
              queryPrimers({ Institution: 'US CDC', Name: '2019-nCoV-N2-F' })
            )
            .concat(
              queryPrimers({ Institution: 'US CDC', Name: '2019-nCoV-N2-R' })
            )
            .concat(
              queryPrimers({ Institution: 'US CDC', Name: '2019-nCoV-N2-P' })
            ),
          startDate: MIN_DATE,
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
          endDate: todayISO(),
        },
      },
    },
    {
      title: 'Co-occurrence of R203K and G204R in N gene',
      description:
<<<<<<< HEAD
        'Two MUTATIONs in the N gene, R203K and G204R, co-occur with each other. These MUTATIONs are associated with the B.1.1 lineage.',
      image: N203204Image,
=======
        'Two mutations in the N gene, R203K and G204R, co-occur with each other. These mutations are associated with the B.1.1 lineage.',
      image:
        'https://storage.googleapis.com/ve-public/example/n_203_204_coocurrence.png',
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
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
<<<<<<< HEAD
          startDate: min_date,
=======
          startDate: MIN_DATE,
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
          endDate: todayISO(),
        },
      },
    },
    {
      title: 'Emergence of Spike D614G in Europe/North America',
      description: '',
<<<<<<< HEAD
      image: D614GEuropeNAImage,
=======
      image:
        'https://storage.googleapis.com/ve-public/example/d614g_europe_na.png',
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
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
<<<<<<< HEAD
          startDate: min_date,
=======
          startDate: MIN_DATE,
>>>>>>> e6dd8312 (Rsvg workflow main (#420))
          endDate: '2020-12-31',
        },
      },
    },
  ];
<<<<<<< HEAD
}
=======
};
>>>>>>> e6dd8312 (Rsvg workflow main (#420))

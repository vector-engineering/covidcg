import React from 'react';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';
import _ from 'underscore';
import useDimensions from 'react-use-dimensions';

import { getGene } from '../../utils/gene_protein';
import { getLocationByNameAndLevel } from '../../utils/location';
import { queryPrimers } from '../../utils/primer';

import {
  GROUP_SNV,
  DNA_OR_AA,
  COORDINATE_MODES,
  NORM_MODES,
  COUNT_MODES,
  DATE_BINS,
  ASYNC_STATES,
  TABS,
} from '../../constants/defs.json';
import { config } from '../../config';

import ExternalLink from '../Common/ExternalLink';
import SkeletonElement from '../Common/SkeletonElement';
// import LoadingSpinner from '../Common/LoadingSpinner';

import SurveillancePlot from '../Vega/SurveillancePlot';
import GlobalSeqPlot from '../Vega/GlobalSeqPlot';

import TempImage from '../../assets/images/cg_short_v13@4x_square.png';

import GlobalLineagesImage from '../../assets/analysis_screens/global_lineages.png';
// import XinfadiImage from '../../assets/analysis_screens/xinfadi.png';
import IcelandImage from '../../assets/analysis_screens/iceland.png';
import D614GWestCoastImage from '../../assets/analysis_screens/d614g_west_coast.png';
import D614GUSStatesImage from '../../assets/analysis_screens/d614g_us_states.png';
import USCDCPrimerImage from '../../assets/analysis_screens/us_cdc_primer.png';
import D614GEuropeNAImage from '../../assets/analysis_screens/d614g_europe_na.png';
import N203204Image from '../../assets/analysis_screens/n_203_204_coocurrence.png';

const ExampleTabContainer = styled.div`
  background-color: #f8f8f8;
`;

const ExampleTabContent = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  flex-grow: 1;
  padding: 20px;
  margin: 0 auto;

  background-color: #fff;
`;

const ExampleHeader = styled.div`
  padding-left: 10px;
  max-width: 800px;

  p {
    font-weight: normal;
    line-height: normal;
  }
`;

const ExampleTitle = styled.h2``;

const ExampleList = styled.div`
  display: flex;
  flex-direction: row;
  flex-wrap: wrap;
`;
const ExampleItem = styled.a`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  width: 300px;
  height: 250px;
  margin: 10px;
  border: 1px solid #ccc;
  box-shadow: 0px 3px 4px #aaa;

  text-decoration: none;

  transition: 0.1s all ease-in-out;
  &:hover {
    border: 1px solid #00e;
  }
`;
const ExampleItemImage = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: center;
  overflow: hidden;

  img {
    width: auto;
    height: 150px;
  }
`;
const ExampleItemFooter = styled.div`
  height: 100px;
  padding: 10px;

  border-top: 1px solid #aaa;

  .example-item-title {
    display: block;
    font-size: 1.1em;
    font-weight: 700;
    color: #444;
  }

  .example-item-description {
    display: block;
    font-size: 0.85em;
    font-weight: normal;
    color: #888;
    margin: 3px 0px;
    line-height: normal;
  }
  //
`;

// const TOC = styled.div`
//   font-weight: normal;
// `;

// const ExampleTutorial = styled.div`
//   padding: 10px;
//   padding-top: 0px;
//   font-weight: normal;
//   line-height: normal;
//   max-width: 800px;
// `;

// const TutorialImage = styled.img`
//   width: 750px;
//   border: 1px solid #aaa;
// `;

// const scrollToRef = (id, e) => {
//   if (e !== undefined) {
//     e.preventDefault();
//   }
//   const el = window.document.getElementById(id);
//   window.scrollTo(0, el.offsetTop);
// };

const ExampleTab = observer(() => {
  const {
    configStore,
    plotSettingsStore,
    UIStore,
    locationDataStore,
  } = useStores();

  const [ref, { width }] = useDimensions();

  const selectTree = locationDataStore.selectTree;
  const exampleItems = [
    {
      title: 'Global Lineages',
      description: 'View the growth of the B lineage family over all locations',
      image: GlobalLineagesImage,
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_PERCENTAGES,
          groupStackCountMode: COUNT_MODES.COUNT_NEW,
          groupStackDateBin: DATE_BINS.DATE_BIN_WEEK,
        },
        UI: {
          activeTab: TABS.TAB_GROUP,
        },
        config: {
          groupKey: config.group_cols.lineage.name,
          dnaOrAa: DNA_OR_AA.DNA,
          selectedGene: getGene('S'),
          coordinateMode: COORDINATE_MODES.COORD_GENE,
          selectedLocationNodes: [selectTree], // select root
        },
      },
    },
    // {
    //   title: 'Lineages in China - Beijing Xinfadi Market',
    //   description:
    //     "New lineages in uncovered in early June in Beijing's Xinfadi Market may have been circulating in China in March",
    //   image: XinfadiImage,
    //   settings: {
    //     plotSettings: {
    //       groupStackNormMode: NORM_MODES.NORM_COUNTS,
    //       groupStackCountMode: COUNT_MODES.COUNT_NEW,
    //       groupStackDateBin: DATE_BINS.DATE_BIN_DAY,
    //     },
    //     UI: {
    //       activeTab: TABS.TAB_GROUP,
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
      image: IcelandImage,
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_COUNTS,
          groupStackCountMode: COUNT_MODES.COUNT_NEW,
          groupStackDateBin: DATE_BINS.DATE_BIN_DAY,
        },
        UI: {
          activeTab: TABS.TAB_GROUP,
        },
        config: {
          groupKey: config.group_cols.lineage.name,
          dnaOrAa: DNA_OR_AA.AA,
          selectedGene: getGene('S'),
          coordinateMode: COORDINATE_MODES.COORD_GENE,
          selectedLocationNodes: [
            getLocationByNameAndLevel(selectTree, 'Iceland', 'country')[0],
          ],
        },
      },
    },
    {
      title: 'Rise of Spike D614G mutation in West Coast USA',
      description:
        'The Spike D614G mutation has accumulated in frequency in the West Coast states of the USA',
      image: D614GWestCoastImage,
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_PERCENTAGES,
          groupStackCountMode: COUNT_MODES.COUNT_NEW,
          groupStackDateBin: DATE_BINS.DATE_BIN_DAY,
        },
        UI: {
          activeTab: TABS.TAB_GROUP,
        },
        config: {
          groupKey: GROUP_SNV,
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
      image: D614GUSStatesImage,
      settings: {
        plotSettings: {
          locationDateNormMode: NORM_MODES.NORM_PERCENTAGES,
          locationDateCountMode: COUNT_MODES.COUNT_CUMULATIVE,
          locationDateDateBin: DATE_BINS.DATE_BIN_MONTH,
        },
        UI: {
          activeTab: TABS.TAB_LOCATION,
        },
        config: {
          groupKey: GROUP_SNV,
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
        },
      },
    },
    {
      title: 'Diagnostics: NT SNVs in US CDC qPCR primer/probe sequences',
      description:
        'Prevalence of any mutations present within the US CDC primer and probe sequences (N1 + N2), for sequences in the US',
      image: USCDCPrimerImage,
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_PERCENTAGES,
          groupStackCountMode: COUNT_MODES.COUNT_CUMULATIVE,
          groupStackDateBin: DATE_BINS.DATE_BIN_DAY,
        },
        UI: {
          activeTab: TABS.TAB_GROUP,
        },
        config: {
          groupKey: GROUP_SNV,
          dnaOrAa: DNA_OR_AA.DNA,
          coordinateMode: COORDINATE_MODES.COORD_PRIMER,
          selectedLocationNodes: [
            getLocationByNameAndLevel(selectTree, 'USA', 'country')[0],
          ],
          selectedPrimers: []
            .concat(
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
        },
      },
    },
    {
      title: 'Co-occurrence of R203K and G204R in N gene',
      description:
        'Two SNVs in the N gene, R203K and G204R, co-occur with each other. These SNVs are associated with the B.1.1 lineage.',
      image: N203204Image,
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_PERCENTAGES,
          groupStackCountMode: COUNT_MODES.COUNT_NEW,
          groupStackDateBin: DATE_BINS.DATE_BIN_WEEK,
        },
        UI: {
          activeTab: TABS.TAB_GROUP,
        },
        config: {
          groupKey: GROUP_SNV,
          dnaOrAa: DNA_OR_AA.AA,
          selectedGene: getGene('N'),
          coordinateMode: COORDINATE_MODES.COORD_GENE,
          selectedLocationNodes: [
            getLocationByNameAndLevel(selectTree, 'USA', 'country')[0],
          ],
          selectedGroups: [{ group: 'N|203|R|K' }, { group: 'N|204|G|R' }],
        },
      },
    },
    {
      title: 'Emergence of Spike D614G in Europe/North America',
      description: '',
      image: D614GEuropeNAImage,
      settings: {
        plotSettings: {
          groupStackNormMode: NORM_MODES.NORM_PERCENTAGES,
          groupStackCountMode: COUNT_MODES.COUNT_NEW,
          groupStackDateBin: DATE_BINS.DATE_BIN_WEEK,
        },
        UI: {
          activeTab: TABS.TAB_GROUP,
        },
        config: {
          groupKey: GROUP_SNV,
          dnaOrAa: DNA_OR_AA.AA,
          selectedGene: getGene('S'),
          coordinateMode: COORDINATE_MODES.COORD_GENE,
          selectedLocationNodes: [
            getLocationByNameAndLevel(selectTree, 'North America', 'region')[0],
            getLocationByNameAndLevel(selectTree, 'Europe', 'region')[0],
          ],
          selectedGroups: [{ group: 'S|614|D|G' }],
        },
      },
    },
  ];

  const onExampleClick = (title, e) => {
    e.preventDefault();

    const exampleItem = _.findWhere(exampleItems, { title });
    // console.log(exampleItem);

    // Apply settings for each store
    Object.keys(exampleItem.settings).forEach((store) => {
      // Call the example action for each store
      if (store === 'config') {
        configStore.resetValues(exampleItem.settings[store]);
      } else if (store === 'plotSettings') {
        plotSettingsStore.resetValues(exampleItem.settings[store]);
      } else if (store === 'UI') {
        UIStore.resetValues(exampleItem.settings[store]);
      }
    });
  };

  const exampleElements = [];
  exampleItems.forEach((exampleItem) => {
    exampleElements.push(
      <ExampleItem
        key={`example-item-${exampleItem.title}`}
        href="#"
        onClick={onExampleClick.bind(this, exampleItem.title)}
      >
        <ExampleItemImage>
          <img
            src={
              exampleItem.image === undefined ? TempImage : exampleItem.image
            }
          />
        </ExampleItemImage>
        <ExampleItemFooter>
          <span className="example-item-title">{exampleItem.title}</span>
          <p className="example-item-description">{exampleItem.description}</p>
        </ExampleItemFooter>
      </ExampleItem>
    );
  });

  const renderExamples = () => {
    // Hide the examples while the app is still initializing
    if (UIStore.caseDataState === ASYNC_STATES.STARTED) {
      const skeletonList = [];
      for (let i = 0; i < exampleElements.length; i++) {
        skeletonList.push(
          <SkeletonElement
            key={`example-loading-${i}`}
            delay={2}
            width="300px"
            height={250}
            style={{
              margin: '10px',
            }}
          />
        );
      }
      return skeletonList;
    }

    return exampleElements;
  };

  return (
    <ExampleTabContainer ref={ref} r>
      <ExampleTabContent>
        {/* <TOC>
          <h3>Table of Contents</h3>
          <ul>
            <li>
              <a
                href="#introduction"
                onClick={scrollToRef.bind(this, 'introduction')}
              >
                <b>Introduction</b>
              </a>
            </li>
            <li>
              <a
                href="#example-analyses"
                onClick={scrollToRef.bind(this, 'example-analyses')}
              >
                <b>Example Analyses</b>
              </a>
            </li>
          </ul>
        </TOC> */}

        <SurveillancePlot width={width - 150} />
        <div style={{ height: '15px' }} />
        <GlobalSeqPlot width={width - 120} />

        <ExampleHeader>
          {/* <a id="example-analyses" /> */}
          <ExampleTitle>Example Analyses</ExampleTitle>
          {/* <a href="#" onClick={scrollToRef.bind(this, 'getting-started-top')}>
            Back
          </a> */}
          <p>
            Use these example analyses to get started and explore the features
            of this application. If you would like to add an analysis to this
            list, please{' '}
            <ExternalLink href="https://github.com/vector-engineering/covidcg">
              submit a pull request on our GitHub
            </ExternalLink>
            , or contact us at{' '}
            <ExternalLink href="mailto:covidcg@broadinstitute.org">
              covidcg@broadinstitute.org
            </ExternalLink>
          </p>
        </ExampleHeader>
        <ExampleList>{renderExamples()}</ExampleList>
      </ExampleTabContent>
    </ExampleTabContainer>
  );
});

export default ExampleTab;

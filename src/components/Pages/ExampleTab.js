import React from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import styled from 'styled-components';
import { useStores } from '../../stores/connect';
import _ from 'underscore';

import { ASYNC_STATES, TABS } from '../../constants/UI';

import { getGene } from '../../utils/gene';
import { getProtein } from '../../utils/protein';
import {
  getLocationIds,
  getLocationByNameAndLevel,
  loadSelectTree,
} from '../../utils/location';
import { queryPrimers } from '../../utils/primer';

import {
  GROUP_KEYS,
  DNA_OR_AA,
  COORDINATE_MODES,
} from '../../constants/config';
import {
  NORM_MODES,
  COUNT_MODES,
  DATE_BINS,
} from '../../constants/plotSettings';

import SkeletonElement from '../Common/SkeletonElement';
import LoadingSpinner from '../Common/LoadingSpinner';
import TempImage from '../../assets/images/cg_short_v13@4x_square.png';

const ExampleTabContainer = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  padding: 20px;
`;

const ExampleHeader = styled.div`
  padding-left: 10px;
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
`;

const selectTree = loadSelectTree();
const NorthAmericaNode = getLocationByNameAndLevel(
  selectTree,
  'North America',
  'region'
)[0];
const EuropeNode = getLocationByNameAndLevel(selectTree, 'Europe', 'region')[0];
const ChinaNode = getLocationByNameAndLevel(selectTree, 'China', 'country')[0];

const WestCoastNodes = [
  getLocationByNameAndLevel(selectTree, 'Washington', 'division')[0],
  getLocationByNameAndLevel(selectTree, 'Oregon', 'division')[0],
  getLocationByNameAndLevel(selectTree, 'California', 'division')[0],
];

const USStateNodes = [
  getLocationByNameAndLevel(selectTree, 'Washington', 'division')[0],
  getLocationByNameAndLevel(selectTree, 'California', 'division')[0],
  getLocationByNameAndLevel(selectTree, 'Florida', 'division')[0],
  getLocationByNameAndLevel(selectTree, 'Massachusetts', 'division')[0],
  getLocationByNameAndLevel(selectTree, 'Michigan', 'division')[0],
  getLocationByNameAndLevel(selectTree, 'New York', 'division')[0],
  getLocationByNameAndLevel(selectTree, 'Texas', 'division')[0],
  getLocationByNameAndLevel(selectTree, 'Wisconsin', 'division')[0],
];

const USANode = getLocationByNameAndLevel(selectTree, 'USA', 'country')[0];

const CDCPrimers = []
  .concat(queryPrimers({ Institution: 'US CDC', Name: '2019-nCoV-N1-F' }))
  .concat(queryPrimers({ Institution: 'US CDC', Name: '2019-nCoV-N1-R' }))
  .concat(queryPrimers({ Institution: 'US CDC', Name: '2019-nCoV-N1-P' }))
  .concat(queryPrimers({ Institution: 'US CDC', Name: '2019-nCoV-N2-F' }))
  .concat(queryPrimers({ Institution: 'US CDC', Name: '2019-nCoV-N2-R' }))
  .concat(queryPrimers({ Institution: 'US CDC', Name: '2019-nCoV-N2-P' }));

const exampleItems = [
  {
    title: 'Global Lineages',
    description: 'View the growth of the B lineage family over all locations',
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
        groupKey: GROUP_KEYS.GROUP_LINEAGE,
        dnaOrAa: DNA_OR_AA.DNA,
        selectedGene: getGene('S'),
        coordinateMode: COORDINATE_MODES.COORD_GENE,
        coordinateRanges: getGene('S').ranges,
        selectedLocationNodes: [selectTree], // select root
        selectedLocationIds: getLocationIds([selectTree]),
      },
    },
  },
  {
    title: 'Lineages in China - Beijing Xinfadi Market',
    description:
      "New lineages in uncovered in early June in Beijing's Xinfadi Market may have been circulating in China in March",
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
        groupKey: GROUP_KEYS.GROUP_LINEAGE,
        dnaOrAa: DNA_OR_AA.AA,
        selectedGene: getGene('S'),
        coordinateMode: COORDINATE_MODES.COORD_GENE,
        coordinateRanges: getGene('S').ranges,
        selectedLocationNodes: [ChinaNode],
        selectedLocationIds: getLocationIds([ChinaNode]),
      },
    },
  },
  {
    title: 'Rise of Spike D614G mutation in West Coast USA',
    description:
      'The Spike D614G mutation has accumulated in frequency in the West Coast states of the USA',
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
        groupKey: GROUP_KEYS.GROUP_SNV,
        dnaOrAa: DNA_OR_AA.AA,
        selectedGene: getGene('S'),
        coordinateMode: COORDINATE_MODES.COORD_GENE,
        coordinateRanges: getGene('S').ranges,
        selectedLocationNodes: WestCoastNodes,
        selectedLocationIds: getLocationIds(WestCoastNodes),
      },
    },
  },
  {
    title: 'Prevalence of Spike D614G in various US States',
    description:
      'The proportion of sequences with the Spike D614G mutation varies between US States',
    settings: {
      plotSettings: {
        locationDateNormMode: NORM_MODES.NORM_PERCENTAGES,
        locationDateCountMode: COUNT_MODES.COUNT_CUMULATIVE,
        locationDateDateBin: DATE_BINS.DATE_BIN_DAY,
      },
      UI: {
        activeTab: TABS.TAB_LOCATION,
      },
      config: {
        groupKey: GROUP_KEYS.GROUP_SNV,
        dnaOrAa: DNA_OR_AA.AA,
        selectedGene: getGene('S'),
        coordinateMode: COORDINATE_MODES.COORD_GENE,
        coordinateRanges: getGene('S').ranges,
        selectedLocationNodes: USStateNodes,
        selectedLocationIds: getLocationIds(USStateNodes),
        selectedGroups: [{ group: 'S|614|D|G' }],
      },
    },
  },
  {
    title: 'Diagnostics: NT SNVs in US CDC qPCR primer/probe sequences',
    description:
      'Prevalence of any mutations present within the US CDC primer and probe sequences (N1 + N2), for sequences in the US',
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
        groupKey: GROUP_KEYS.GROUP_SNV,
        dnaOrAa: DNA_OR_AA.DNA,
        coordinateMode: COORDINATE_MODES.COORD_PRIMER,
        coordinateRanges: CDCPrimers.map((primer) => [
          primer.Start,
          primer.End,
        ]),
        selectedLocationNodes: [USANode],
        selectedLocationIds: getLocationIds([USANode]),
        selectedPrimers: CDCPrimers,
      },
    },
  },
  {
    title: 'Emergence of Spike D614G in Europe/North America',
    description: '',
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
        groupKey: GROUP_KEYS.GROUP_SNV,
        dnaOrAa: DNA_OR_AA.AA,
        selectedGene: getGene('S'),
        coordinateMode: COORDINATE_MODES.COORD_GENE,
        coordinateRanges: getGene('S').ranges,
        selectedLocationNodes: [NorthAmericaNode, EuropeNode],
        selectedLocationIds: getLocationIds([NorthAmericaNode, EuropeNode]),
      },
    },
  },
];

const ExampleTab = observer(() => {
  const { configStore, plotSettingsStore, UIStore } = useStores();

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
          <img src={TempImage} />
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
    <ExampleTabContainer>
      <ExampleHeader>
        <ExampleTitle>Example Analyses</ExampleTitle>
        <p>Get started with a bunch of example analyses</p>
      </ExampleHeader>
      <ExampleList>{renderExamples()}</ExampleList>
    </ExampleTabContainer>
  );
});
ExampleTab.propTypes = {
  width: PropTypes.number,
};
ExampleTab.defaultProps = {
  width: 100,
};

export default ExampleTab;

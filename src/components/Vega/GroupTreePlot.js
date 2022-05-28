import React, { useState, useEffect, useRef, useLayoutEffect } from 'react';
import PropTypes from 'prop-types';
import { config } from '../../config';

import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { ISOToInt } from '../../utils/date';
import { TREE_COLOR_MODES } from '../../constants/defs.json';

import QuestionButton from '../Buttons/QuestionButton';
import DropdownButton from '../Buttons/DropdownButton';
import VegaEmbed from '../../react_vega/VegaEmbed';
import ExternalLink from '../Common/ExternalLink';
import GradientLegend from './GradientLegend';
import MarkLegend from './MarkLegend';

import {
  TreePlotContainer,
  Header,
  HeaderRow,
  LegendRow,
  Title,
  SubTitle,
  SelectContainer,
  TreeScrollContainer,
} from './GroupTreePlot.styles';

// import legendSpec from '../../vega_specs/group_tree_legend_v2.vg.json';
import treeSpec from '../../vega_specs/group_tree_v2.vg.json';

const headerHeight = 135;
const treePlotHeight = 12000;

// https://cssgradient.io/
const VIRIDIS_GRADIENT =
  'linear-gradient(90deg, rgba(68,1,84,1) 0%, rgba(57,86,140,1) 25%, rgba(31,150,139,1) 50%, rgba(115,208,85,1) 75%, rgba(253,231,37,1) 100%);';
const INFERNO_GRADIENT =
  'linear-gradient(90deg, rgba(0,0,0,1) 0%, rgba(26,11,64,1) 11%, rgba(74,11,106,1) 22%, rgba(120,28,109,1) 33%, rgba(164,44,96,1) 44%, rgba(207,68,70,1) 55%, rgba(237,104,37,1) 66%, rgba(251,155,6,1) 77%, rgba(247,209,60,1) 88%, rgba(252,254,164,1) 100%);';

const DOWNLOAD_OPTIONS = {
  JSON: 'JSON',
  NEWICK: 'Newick',
};

const GroupTreePlot = observer(({ width }) => {
  // const vegaLegendRef = useRef();
  const vegaRef = useRef();
  const treeContainerRef = useRef();
  const hiddenLink = useRef();
  const [treeContainerDimensions, setTreeContainerDimensions] = useState();
  const { groupDataStore, plotSettingsStore } = useStores();

  // Pulled from react-use-dimensions
  // https://swizec.com/blog/usedimensions-a-react-hook-to-measure-dom-nodes/
  // For some reason, the ref packaged with the
  // useDimensions() hook does not work with scrolling as well
  useLayoutEffect(() => {
    setTreeContainerDimensions(
      treeContainerRef.current.getBoundingClientRect().toJSON()
    );
  }, [treeContainerRef.current]);

  const handleSelected = (...args) => {
    // console.log(args);

    // Don't fire if the selection is the same
    if (
      args[1].length === groupDataStore.selectedGroups.length &&
      args[1]
        .map((group) => group.lineage)
        .every((group) => groupDataStore.selectedGroups.includes(group))
    ) {
      return;
    }

    groupDataStore.updateSelectedGroups(args[1].map((group) => group.lineage));
  };

  const processSelectedGroups = () => {
    return JSON.parse(
      JSON.stringify(
        groupDataStore.selectedGroups.map((group) => {
          return { lineage: group };
        })
      )
    );
  };

  const [state, setState] = useState({
    data: {
      selected: processSelectedGroups(),
    },
    dataListeners: {
      selected: handleSelected,
    },
  });

  const onChangeTreeColorMode = (event) => {
    plotSettingsStore.setReportTreeColorMode(event.target.value);
  };

  const handleDownloadSelect = (option) => {
    if (option === DOWNLOAD_OPTIONS.JSON) {
      hiddenLink.current.href =
        'https://storage.googleapis.com/ve-public/lineage_tree_graph.json';
      hiddenLink.current.click();
    } else if (option === DOWNLOAD_OPTIONS.NEWICK) {
      hiddenLink.current.href =
        'https://storage.googleapis.com/ve-public/lineage_representatives_cleaned.nwk';
      hiddenLink.current.click();
    }
  };

  // Update internal selected groups copy
  useEffect(() => {
    vegaRef.current.getData('selected', (selected) => {
      // If the internal Vega copy of the selected groups is smaller
      // than the store's version, that means we just selected a new group
      // from outside of the Vega plot
      // Use this to trigger a zoom to the newly selected group
      if (selected.length < groupDataStore.selectedGroups.length) {
        const newGroup = groupDataStore.selectedGroups.find(
          (group) => !selected.map((g) => g.lineage).includes(group)
        );
        // console.log(newGroup);

        vegaRef.current.getData('tree', (tree) => {
          // console.log(treeContainerDimensions);
          const treeContainerHeight = treeContainerDimensions.height;
          // Find the y-position of the newly selected lineage
          let scrollToY = tree.find((row) => row.lineage === newGroup).y;
          // Subtract half of the height of the tree container so that it's positioned
          // in the middle of the plot
          scrollToY -= treeContainerHeight / 2;
          // Enforce upper bound
          if (scrollToY < 0) {
            scrollToY = 0;
          }
          // Enforce lower bound
          // Not sure if this is necessary for all browsers?
          else if (scrollToY > treePlotHeight - treeContainerHeight) {
            scrollToY = treePlotHeight - treeContainerHeight;
          }

          treeContainerRef.current.scroll({
            top: scrollToY,
            left: 0,
            behavior: 'smooth', // or 'auto'
          });
        });
      }
    });

    setState({
      ...state,
      data: {
        ...state.data,
        selected: processSelectedGroups(),
      },
    });
  }, [groupDataStore.selectedGroups]);

  const renderLegend = () => {
    if (
      plotSettingsStore.reportTreeColorMode === TREE_COLOR_MODES.COLOR_REGION
    ) {
      return (
        <MarkLegend
          title={'Most Common Region'}
          itemLabels={[
            'Africa',
            'Asia',
            'Europe',
            'North America',
            'Oceania',
            'South America',
          ]}
          itemColors={[
            // TABLEAU10
            '#4c78a8',
            '#f58518',
            '#e45756',
            '#72b7b2',
            '#54a24b',
            '#eeca3b',
            '#b279a2',
            '#ff9da6',
            '#9d755d',
            '#bab0ac',
          ]}
        />
      );
    } else if (
      plotSettingsStore.reportTreeColorMode === TREE_COLOR_MODES.COLOR_LATEST
    ) {
      // A tick for each 6 months?
      const startDate = ISOToInt(config.min_date);
      const endDate = Date.now();
      // Overall date range in milliseconds
      const rangeMS = endDate - startDate;

      // Convert ISO strings to fractional tick locations
      let ticks = ['2020-01-01', '2020-07-01', '2021-01-01', '2021-07-01'];
      let tickLocs = ticks.map((tick) => {
        return (ISOToInt(tick) - startDate) / rangeMS;
      });
      let tickLabels = ['2020-01', '07', '2021-01', '07'];

      return (
        <GradientLegend
          title={'Most Recent Collection Date'}
          gradient={VIRIDIS_GRADIENT}
          ticks={tickLocs}
          tickLabels={tickLabels}
        />
      );
    } else if (
      plotSettingsStore.reportTreeColorMode ===
      TREE_COLOR_MODES.COLOR_NUM_SEQUENCES
    ) {
      let tickLocs = [1, 10, 100, 1000, 10000, 100000];
      let tickLabels = tickLocs.map((tick) => tick.toString());
      tickLocs = tickLocs.map((tick) => Math.log10(tick) / Math.log10(100000));

      return (
        <GradientLegend
          title={'Num. Sequences'}
          gradient={INFERNO_GRADIENT}
          ticks={tickLocs}
          tickLabels={tickLabels}
        />
      );
    }
  };

  return (
    <TreePlotContainer width={width}>
      <a
        ref={hiddenLink}
        href=""
        target="_blank"
        rel="noopener noreferrer"
        style={{ visibility: 'hidden' }}
      />
      <Header headerHeight={headerHeight}>
        <HeaderRow>
          <Title>Time-Scaled Tree</Title>
          <QuestionButton
            data-tip='<p>Displayed below is a phylogenetic tree of PANGO lineages (each lineage represented by a "representative genome".</p><p>Branch points in the tree indicate the <i>predicted</i> date of divergence.</p><p>The bars represent the collection range for each lineage. The start of the bar represents the earliest collection date, and the end of the bar represents the most recent collection date.</p><p>In the dropdown below, you can color the bars by collection date, # sequences, or by the region (continent) with the most counts of that lineage.</p>'
            data-html="true"
            data-for="main-tooltip"
          />
        </HeaderRow>
        <HeaderRow>
          <SubTitle>
            Methodology and Visualization from{' '}
            <ExternalLink href="https://filogeneti.ca/CoVizu/">
              CoVizu
            </ExternalLink>
          </SubTitle>
        </HeaderRow>
        <HeaderRow>
          <DropdownButton
            text={'Download'}
            options={[DOWNLOAD_OPTIONS.JSON, DOWNLOAD_OPTIONS.NEWICK]}
            onSelect={handleDownloadSelect}
            style={{
              padding: '2px 10px',
              paddingRight: '18px',
              marginBottom: '2px',
            }}
            direction={'left'}
          />
        </HeaderRow>
        <SelectContainer>
          <label>
            Color tree by
            <select
              value={plotSettingsStore.reportTreeColorMode}
              onChange={onChangeTreeColorMode}
            >
              <option value={TREE_COLOR_MODES.COLOR_REGION}>Region</option>
              <option value={TREE_COLOR_MODES.COLOR_LATEST}>
                Collection Date
              </option>
              <option value={TREE_COLOR_MODES.COLOR_NUM_SEQUENCES}>
                # Sequences
              </option>
            </select>
          </label>
        </SelectContainer>
        <LegendRow>{renderLegend()}</LegendRow>
      </Header>
      {/* <VegaEmbed
        ref={vegaLegendRef}
        spec={legendSpec}
        width={width - 75}
        actions={false}
      /> */}
      <TreeScrollContainer ref={treeContainerRef} headerHeight={headerHeight}>
        <VegaEmbed
          ref={vegaRef}
          data={state.data}
          spec={treeSpec}
          // signalListeners={state.signalListeners}
          dataListeners={state.dataListeners}
          signals={{
            plotHeight: treePlotHeight,
            colorScale: plotSettingsStore.reportTreeColorMode,
            latestDomain: [ISOToInt('2019-12-15'), Date.now()],
            numSequencesDomain: [1, 100000],
          }}
          width={width - 75}
          actions={false}
        />
      </TreeScrollContainer>
    </TreePlotContainer>
  );
});

GroupTreePlot.propTypes = {
  width: PropTypes.number.isRequired,
};
GroupTreePlot.defaultProps = {
  width: 400,
};

export default GroupTreePlot;

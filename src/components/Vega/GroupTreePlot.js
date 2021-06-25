import React, { useState, useEffect, useRef, useLayoutEffect } from 'react';
import PropTypes from 'prop-types';

import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { TREE_COLOR_MODES } from '../../constants/defs.json';

import VegaEmbed from '../../react_vega/VegaEmbed';
import ExternalLink from '../Common/ExternalLink';

import {
  TreePlotContainer,
  Header,
  Title,
  SubTitle,
  SelectContainer,
  TreeScrollContainer,
} from './GroupTreePlot.styles';

// import legendSpec from '../../vega_specs/group_tree_legend_v2.vg.json';
import treeSpec from '../../vega_specs/group_tree_v2.vg.json';

const headerHeight = 60;
const treePlotHeight = 12000;

const GroupTreePlot = observer(({ width }) => {
  // const vegaLegendRef = useRef();
  const vegaRef = useRef();
  const treeContainerRef = useRef();
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

  return (
    <TreePlotContainer width={width}>
      <Header headerHeight={headerHeight}>
        <Title>Time-Scaled Tree</Title>
        <SubTitle>
          Methodology and Visualization from{' '}
          <ExternalLink href="https://filogeneti.ca/CoVizu/">
            CoVizu
          </ExternalLink>
        </SubTitle>
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

import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import VegaEmbed from '../../react_vega/VegaEmbed';
import initialSpec from '../../vega_specs/group_tree.vg.json';

const headerHeight = 40;

const TreePlotContainer = styled.div`
  width: ${({ width }) => width}px;
  height: 100vh;

  display: flex;
  flex-direction: column;
  align-items: stretch;

  position: sticky;
  left: 0px;
  top: 0px;
`;
TreePlotContainer.defaultProps = {
  width: 300,
};

const Header = styled.div`
  height: ${({ headerHeight }) => headerHeight}px;
`;

const TreeScrollContainer = styled.div`
  overflow-y: scroll;
  height: calc(100vh - ${({ headerHeight }) => headerHeight}px);
`;

const GroupTreePlot = observer(({ width }) => {
  const vegaRef = useRef();
  // const { dataStore, UIStore, configStore, plotSettingsStore } = useStores();

  return (
    <TreePlotContainer width={width}>
      <Header headerHeight={headerHeight}>
        <span>Phylogenetic Tree</span>
      </Header>
      <TreeScrollContainer headerHeight={headerHeight}>
        <VegaEmbed
          ref={vegaRef}
          // data={state.data}
          spec={initialSpec}
          // signalListeners={state.signalListeners}
          // dataListeners={state.dataListeners}
          signals={
            {
              // disableSelectionColoring: configStore.groupKey === GROUP_SNV,
              // detailHeight,
              // hoverBar: { group: configStore.hoverGroup },
              // stackOffset,
              // dateBin,
              // cumulativeWindow,
              // detailYLabel,
              // yFormat:
              //   plotSettingsStore.groupStackNormMode === NORM_MODES.NORM_COUNTS
              //     ? 's'
              //     : '%',
            }
          }
          cheapSignals={['hoverBar']}
          width={width - 75}
          height={2000}
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

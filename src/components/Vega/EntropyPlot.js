import React, { useState, useEffect } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { toJS } from 'mobx';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import _ from 'underscore';

import VegaEmbed from '../../react_vega/VegaEmbed';

import initialSpec from '../../vega_specs/entropy.vg.json';

const PlotContainer = styled.div``;

const EntropyPlot = observer(({ width }) => {
  const { configStore, dataStore } = useStores();

  const processData = (groupCountArr) => {
    return groupCountArr;
  };

  const handleHoverGroup = (...args) => {
    // Don't fire the action if there's no change
    let hoverGroup = args[1] === null ? null : args[1]['group'];
    if (hoverGroup === configStore.hoverGroup) {
      return;
    }
    configStore.updateHoverGroup(hoverGroup);
  };

  const handleSelected = (...args) => {
    // Don't fire if the selection is the same
    if (_.isEqual(args[1], configStore.selectedGroups)) {
      return;
    }
    configStore.updateSelectedGroups(args[1]);
  };

  const [state, setState] = useState({
    data: {
      table: processData(dataStore.groupCountArr),
      selected: JSON.parse(JSON.stringify(configStore.selectedGroups)),
    },
    signalListeners: {
      hoverGroup: _.throttle(handleHoverGroup, 100),
    },
    dataListeners: {
      selected: handleSelected,
    },
  });

  useEffect(() => {
    setState({
      ...state,
      data: {
        ...state.data,
        table: processData(dataStore.groupCountArr),
      },
    });
  }, [dataStore.selectedRowsAndDateHash]);

  // Update internal selected groups copy
  useEffect(() => {
    setState({
      ...state,
      data: {
        ...state.data,
        selected: JSON.parse(JSON.stringify(configStore.selectedGroups)),
      },
    });
  }, [configStore.selectedGroups]);

  return (
    <PlotContainer>
      <VegaEmbed
        spec={initialSpec}
        data={state.data}
        width={width}
        signals={{
          xLabel: configStore.coordinateMode,
          hoverGroup: { group: configStore.hoverGroup },
        }}
        signalListeners={state.signalListeners}
        dataListeners={state.dataListeners}
      />
    </PlotContainer>
  );
});
EntropyPlot.propTypes = {
  width: PropTypes.number,
};
EntropyPlot.defaultProps = {
  width: 100,
};

export default EntropyPlot;

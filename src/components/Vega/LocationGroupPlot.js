import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import _ from 'underscore';

import VegaEmbed from '../../react_vega/VegaEmbed';
import initialSpec from '../../vega_specs/location_group.vg.json';

const PlotContainer = styled.div``;

const LocationGroupPlot = observer(({ width }) => {
  const vegaRef = useRef();
  const { covidStore } = useStores();

  const handleHoverLocation = (...args) => {
    // Don't fire the action if there's no change
    let hoverLocation = args[1] === null ? null : args[1]['location'];
    if (hoverLocation === covidStore.hoverLocation) {
      return;
    }
    covidStore.updateHoverLocation(hoverLocation);
  };

  const handleHoverGroup = (...args) => {
    // Don't fire the action if there's no change
    let hoverGroup = args[1] === null ? null : args[1]['group'];
    if (hoverGroup === covidStore.hoverGroup) {
      return;
    }
    covidStore.updateHoverGroup(hoverGroup);
  };

  const handleSelectedLocations = (...args) => {
    // console.log(args);
    // Don't fire if the selection is the same
    if (_.isEqual(args[1], covidStore.focusedLocations)) {
      return;
    }
    covidStore.updateFocusedLocations(args[1]);
  };

  const handleSelectedGroups = (...args) => {
    // Don't fire if the selection is the same
    if (_.isEqual(args[1], covidStore.selectedGroups)) {
      return;
    }
    covidStore.updateSelectedGroups(args[1]);
  };

  const [state, setState] = useState({
    data: {
      location_by_group: [],
      selectedGroups: [],
      selectedLocations: JSON.parse(
        JSON.stringify(covidStore.focusedLocations)
      ),
    },
    spec: JSON.parse(JSON.stringify(initialSpec)),
    signalListeners: {
      hoverLocation: _.throttle(handleHoverLocation, 100),
      hoverGroup: _.throttle(handleHoverGroup, 100),
    },
    dataListeners: {
      selectedLocations: handleSelectedLocations,
      selectedGroups: handleSelectedGroups,
    },
  });

  useEffect(() => {
    setState({
      ...state,
      data: {
        ...state.data,
        selectedLocations: JSON.parse(
          JSON.stringify(covidStore.focusedLocations)
        ),
      },
    });
  }, [covidStore.focusedLocations]);

  useEffect(() => {
    setState({
      ...state,
      data: {
        ...state.data,
        location_by_group: JSON.parse(
          JSON.stringify(covidStore.aggLocationData)
        ),
        selectedGroups: JSON.parse(JSON.stringify(covidStore.selectedGroups)),
      },
    });
  }, [covidStore.aggLocationData, covidStore.selectedGroups]);

  let xLabel = 'Sequences by ';
  if (covidStore.groupKey === 'lineage') {
    xLabel += 'Lineage ';
  } else if (covidStore.groupKey === 'clade') {
    xLabel += 'Clade ';
  } else if (covidStore.groupKey === 'snp') {
    if (covidStore.dnaOrAa === 'dna') {
      xLabel += 'NT';
    } else {
      xLabel += 'AA';
    }
    xLabel += ' SNV ';
  }
  xLabel += ' (Cumulative, All Sequences)';

  return (
    <PlotContainer>
      <div style={{ width: `${width}` }}>
        <VegaEmbed
          ref={vegaRef}
          data={state.data}
          spec={state.spec}
          signalListeners={state.signalListeners}
          dataListeners={state.dataListeners}
          signals={{
            hoverLocation: { location: covidStore.hoverLocation },
            hoverGroup: { group: covidStore.hoverGroup },
            xLabel,
          }}
          width={width}
          actions={false}
        />
      </div>
    </PlotContainer>
  );
});
LocationGroupPlot.propTypes = {
  width: PropTypes.number,
};
LocationGroupPlot.defaultProps = {
  width: 100,
};

export default LocationGroupPlot;

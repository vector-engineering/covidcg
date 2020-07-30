import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import _ from 'underscore';

import EmptyPlot from '../Common/EmptyPlot';
import VegaEmbed from '../../react_vega/VegaEmbed';
import initialSpec from '../../vega_specs/location_group.vg.json';

const PlotContainer = styled.div``;

const LocationGroupPlot = observer(({ width }) => {
  const vegaRef = useRef();
  const { dataStore } = useStores();

  const handleHoverLocation = (...args) => {
    // Don't fire the action if there's no change
    let hoverLocation = args[1] === null ? null : args[1]['location'];
    if (hoverLocation === dataStore.hoverLocation) {
      return;
    }
    dataStore.updateHoverLocation(hoverLocation);
  };

  const handleHoverGroup = (...args) => {
    // Don't fire the action if there's no change
    let hoverGroup = args[1] === null ? null : args[1]['group'];
    if (hoverGroup === dataStore.hoverGroup) {
      return;
    }
    dataStore.updateHoverGroup(hoverGroup);
  };

  const handleSelectedLocations = (...args) => {
    // console.log(args);
    // Don't fire if the selection is the same
    if (_.isEqual(args[1], dataStore.focusedLocations)) {
      return;
    }
    dataStore.updateFocusedLocations(args[1]);
  };

  const handleSelectedGroups = (...args) => {
    // Don't fire if the selection is the same
    if (_.isEqual(args[1], dataStore.selectedGroups)) {
      return;
    }
    dataStore.updateSelectedGroups(args[1]);
  };

  const [state, setState] = useState({
    data: {
      location_by_group: [],
      selectedGroups: [],
      selectedLocations: JSON.parse(JSON.stringify(dataStore.focusedLocations)),
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
          JSON.stringify(dataStore.focusedLocations)
        ),
      },
    });
  }, [dataStore.focusedLocations]);

  useEffect(() => {
    setState({
      ...state,
      data: {
        ...state.data,
        location_by_group: JSON.parse(
          JSON.stringify(dataStore.aggLocationData)
        ),
        selectedGroups: JSON.parse(JSON.stringify(dataStore.selectedGroups)),
      },
    });
  }, [dataStore.aggLocationData, dataStore.selectedGroups]);

  let xLabel = 'Sequences by ';
  if (dataStore.groupKey === 'lineage') {
    xLabel += 'Lineage ';
  } else if (dataStore.groupKey === 'clade') {
    xLabel += 'Clade ';
  } else if (dataStore.groupKey === 'snp') {
    if (dataStore.dnaOrAa === 'dna') {
      xLabel += 'NT';
    } else {
      xLabel += 'AA';
    }
    xLabel += ' SNV ';
  }
  xLabel += ' (Cumulative, All Sequences)';

  const renderPlot = () => {
    if (dataStore.selectedLocationIds.length == 0) {
      return (
        <EmptyPlot height={100}>
          <p>
            No locations selected. Please select one or more locations from the
            sidebar, under &quot;Selected Locations&quot;, to compare counts of{' '}
            <b>{dataStore.getGroupLabel()}</b> between them.
          </p>
        </EmptyPlot>
      );
    }

    return (
      <div style={{ width: `${width}` }}>
        <VegaEmbed
          ref={vegaRef}
          data={state.data}
          spec={state.spec}
          signalListeners={state.signalListeners}
          dataListeners={state.dataListeners}
          signals={{
            hoverLocation: { location: dataStore.hoverLocation },
            hoverGroup: { group: dataStore.hoverGroup },
            xLabel,
          }}
          width={width}
          actions={false}
        />
      </div>
    );
  };

  return <PlotContainer>{renderPlot()}</PlotContainer>;
});
LocationGroupPlot.propTypes = {
  width: PropTypes.number,
};
LocationGroupPlot.defaultProps = {
  width: 100,
};

export default LocationGroupPlot;

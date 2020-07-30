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
  const { dataStore, configStore } = useStores();

  const handleHoverLocation = (...args) => {
    // Don't fire the action if there's no change
    let hoverLocation = args[1] === null ? null : args[1]['location'];
    if (hoverLocation === configStore.hoverLocation) {
      return;
    }
    configStore.updateHoverLocation(hoverLocation);
  };

  const handleHoverGroup = (...args) => {
    // Don't fire the action if there's no change
    let hoverGroup = args[1] === null ? null : args[1]['group'];
    if (hoverGroup === configStore.hoverGroup) {
      return;
    }
    configStore.updateHoverGroup(hoverGroup);
  };

  const handleSelectedLocations = (...args) => {
    // console.log(args);
    // Don't fire if the selection is the same
    if (_.isEqual(args[1], configStore.focusedLocations)) {
      return;
    }
    configStore.updateFocusedLocations(args[1]);
  };

  const handleSelectedGroups = (...args) => {
    // Don't fire if the selection is the same
    if (_.isEqual(args[1], configStore.selectedGroups)) {
      return;
    }
    configStore.updateSelectedGroups(args[1]);
  };

  const [state, setState] = useState({
    data: {
      location_by_group: [],
      selectedGroups: [],
      selectedLocations: JSON.parse(
        JSON.stringify(configStore.focusedLocations)
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
          JSON.stringify(configStore.focusedLocations)
        ),
      },
    });
  }, [configStore.focusedLocations]);

  useEffect(() => {
    setState({
      ...state,
      data: {
        ...state.data,
        location_by_group: JSON.parse(
          JSON.stringify(dataStore.aggLocationData)
        ),
        selectedGroups: JSON.parse(JSON.stringify(configStore.selectedGroups)),
      },
    });
  }, [dataStore.aggLocationData, configStore.selectedGroups]);

  let xLabel = 'Sequences by ';
  if (configStore.groupKey === 'lineage') {
    xLabel += 'Lineage ';
  } else if (configStore.groupKey === 'clade') {
    xLabel += 'Clade ';
  } else if (configStore.groupKey === 'snp') {
    if (configStore.dnaOrAa === 'dna') {
      xLabel += 'NT';
    } else {
      xLabel += 'AA';
    }
    xLabel += ' SNV ';
  }
  xLabel += ' (Cumulative, All Sequences)';

  const renderPlot = () => {
    if (configStore.selectedLocationIds.length == 0) {
      return (
        <EmptyPlot height={100}>
          <p>
            No locations selected. Please select one or more locations from the
            sidebar, under &quot;Selected Locations&quot;, to compare counts of{' '}
            <b>{configStore.getGroupLabel()}</b> between them.
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
            hoverLocation: { location: configStore.hoverLocation },
            hoverGroup: { group: configStore.hoverGroup },
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

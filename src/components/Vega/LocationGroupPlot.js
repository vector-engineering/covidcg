import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { aggregate } from '../../utils/transform';
import _ from 'underscore';

import EmptyPlot from '../Common/EmptyPlot';
import VegaEmbed from '../../react_vega/VegaEmbed';
import SkeletonElement from '../Common/SkeletonElement';

import { GROUP_KEYS, DNA_OR_AA } from '../../constants/config';
import { GROUPS } from '../../constants/groups';
import { ASYNC_STATES } from '../../constants/UI';
import initialSpec from '../../vega_specs/location_group.vg.json';

const PlotContainer = styled.div``;

const LocationGroupPlot = observer(({ width }) => {
  const vegaRef = useRef();
  const { dataStore, configStore, UIStore } = useStores();

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

  const processLocationByGroup = () => {
    let locationData = JSON.parse(
      JSON.stringify(dataStore.dataAggLocationGroupDate)
    );

    locationData.forEach((row) => {
      if (!dataStore.groupsToKeep.includes(row.group)) {
        row.group = GROUPS.OTHER_GROUP;
        row.groupName = GROUPS.OTHER_GROUP;
        row.color = '#aaa';
      }
    });

    if (configStore.groupKey === GROUP_KEYS.GROUP_SNV) {
      // Filter out 'Reference' group, when in SNV mode
      locationData = locationData.filter((row) => {
        return row.group !== GROUPS.REFERENCE_GROUP;
      });
    }

    // Filter by date
    if (configStore.dateRange[0] != -1 || configStore.dateRange[1] != -1) {
      locationData = locationData.filter((row) => {
        return (
          (configStore.dateRange[0] == -1 ||
            row.date > configStore.dateRange[0]) &&
          (configStore.dateRange[1] == -1 ||
            row.date < configStore.dateRange[1])
        );
      });
    }

    locationData = aggregate({
      data: locationData,
      groupby: ['location', 'date', 'group', 'groupName'],
      fields: ['cases_sum', 'color', 'location_counts'],
      ops: ['sum', 'first', 'max'],
      as: ['cases_sum', 'color', 'location_counts'],
    });

    return locationData;
  };

  const processSelectedGroups = () => {
    return JSON.parse(JSON.stringify(configStore.selectedGroups));
  };

  const processSelectedLocations = () => {
    return JSON.parse(JSON.stringify(configStore.focusedLocations));
  };

  const [state, setState] = useState({
    data: {
      location_by_group: processLocationByGroup(),
      selectedGroups: processSelectedGroups(),
      selectedLocations: processSelectedLocations(),
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
        selectedLocations: processSelectedLocations(),
      },
    });
  }, [configStore.focusedLocations]);

  useEffect(() => {
    if (UIStore.caseDataState !== ASYNC_STATES.SUCCEEDED) {
      return;
    }

    setState({
      ...state,
      data: {
        ...state.data,
        location_by_group: processLocationByGroup(),
        selectedGroups: processSelectedGroups(),
      },
    });
  }, [
    UIStore.caseDataState,
    configStore.selectedGroups,
    dataStore.groupsToKeep,
    configStore.dateRange,
  ]);

  if (UIStore.caseDataState === ASYNC_STATES.STARTED) {
    return (
      <div
        style={{
          paddingTop: '12px',
          paddingRight: '24px',
          paddingLeft: '12px',
          paddingBottom: '24px',
        }}
      >
        <SkeletonElement delay={2} height={100} />
      </div>
    );
  }

  if (configStore.selectedLocationNodes.length == 0) {
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

  let xLabel, xLabelFormat, stackOffset;
  if (configStore.groupKey === GROUP_KEYS.GROUP_LINEAGE) {
    xLabel += 'Lineage ';
  } else if (configStore.groupKey === GROUP_KEYS.GROUP_CLADE) {
    xLabel += 'Clade ';
  } else if (configStore.groupKey === GROUP_KEYS.GROUP_SNV) {
    if (configStore.dnaOrAa === DNA_OR_AA.DNA) {
      xLabel += 'NT';
    } else {
      xLabel += 'AA';
    }
    xLabel += ' SNV ';
  }
  xLabel += ' (Cumulative, All Sequences)';

  if (configStore.groupKey === GROUP_KEYS.GROUP_SNV) {
    xLabelFormat = 's';
    stackOffset = 'zero';
    xLabel = `# Sequences with ${configStore.getGroupLabel()} (Cumulative, All Sequences)`;
  } else {
    xLabelFormat = '%';
    stackOffset = 'normalize';
    xLabel = `% Sequences by ${configStore.getGroupLabel()} (Cumulative, All Sequences)`;
  }

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
            hoverLocation: { location: configStore.hoverLocation },
            hoverGroup: { group: configStore.hoverGroup },
            xLabel,
            xLabelFormat,
            stackOffset,
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

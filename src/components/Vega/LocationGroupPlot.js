import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { toJS } from 'mobx';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { aggregate } from '../../utils/transform';
import { throttle, debounce } from '../../utils/func';

import EmptyPlot from '../Common/EmptyPlot';
import VegaEmbed from '../../react_vega/VegaEmbed';
import SkeletonElement from '../Common/SkeletonElement';
import { PlotOptions, OptionCheckboxContainer } from './Plot.styles';

import { config } from '../../config';
import {
  GROUP_SNV,
  DNA_OR_AA,
  GROUPS,
  ASYNC_STATES,
} from '../../constants/defs.json';
import initialSpec from '../../vega_specs/location_group.vg.json';

const PlotContainer = styled.div``;

const LocationGroupPlot = observer(({ width }) => {
  const vegaRef = useRef();
  const {
    dataStore,
    configStore,
    UIStore,
    plotSettingsStore,
    groupDataStore,
    snpDataStore,
  } = useStores();

  const handleHoverLocation = (...args) => {
    // Don't fire the action if there's no change
    let hoverLocation = args[1] === null ? null : args[1]['location'];
    if (hoverLocation === configStore.hoverLocation) {
      return;
    }
    configStore.updateHoverLocation(hoverLocation);
  };

  const handleHoverGroup = (...args) => {
    configStore.updateHoverGroup(args[1] === null ? null : args[1]['group']);
  };

  const handleSelectedLocations = (...args) => {
    configStore.updateFocusedLocations(args[1]);
  };

  const handleSelectedGroups = (...args) => {
    configStore.updateSelectedGroups(
      args[1] === null
        ? []
        : args[1].map((item) => {
            return { group: item.group };
          })
    );
  };

  const onChangeHideReference = (e) => {
    plotSettingsStore.setLocationGroupHideReference(e.target.checked);
  };

  const processLocationByGroup = () => {
    console.log('LOCATION GROUP PLOT PROCESS DATA');

    let locationData;
    if (configStore.groupKey === GROUP_SNV) {
      locationData = aggregate({
        data: toJS(dataStore.aggLocationSingleSnvDate),
        groupby: ['location', 'group_id'],
        fields: ['counts'],
        ops: ['sum'],
        as: ['counts'],
      });

      if (plotSettingsStore.locationGroupHideReference) {
        // Filter out 'Reference' group, when in SNV mode
        locationData = locationData.filter((row) => {
          return row.group_id !== GROUPS.REFERENCE_GROUP;
        });
      }

      locationData.forEach((record) => {
        let snv = snpDataStore.intToSnv(
          configStore.dnaOrAa,
          configStore.coordinateMode,
          record.group_id
        );
        record.color = snv.color;
        record.group = snv.snp_str;
        record.group_name = snv.name;
      });
    } else {
      locationData = aggregate({
        data: toJS(dataStore.aggLocationGroupDate),
        groupby: ['location', 'group_id'],
        fields: ['counts'],
        ops: ['sum'],
        as: ['counts'],
      });
      locationData.forEach((record) => {
        record.color = groupDataStore.getGroupColor(
          configStore.groupKey,
          record.group_id
        );
        record.group = record.group_id;
        record.group_name = record.group_id;
      });
    }

    locationData.forEach((record) => {
      record.location_counts = dataStore.countsPerLocationMap[record.location];
    });

    // console.log(JSON.stringify(locationData));

    return locationData;
  };

  const [state, setState] = useState({
    // data: {
    //   location_by_group: [],
    //   selectedGroups: [],
    //   selectedLocations: [],
    // },
    hoverGroup: null,
    hoverLocation: null,
    spec: JSON.parse(JSON.stringify(initialSpec)),
    signalListeners: {
      hoverLocation: throttle(handleHoverLocation, 100),
      hoverGroup: debounce(handleHoverGroup, 20),
    },
    dataListeners: {
      selectedLocations: handleSelectedLocations,
      selectedGroups: handleSelectedGroups,
    },
  });

  useEffect(() => {
    setState({
      ...state,
      hoverGroup: { group: configStore.hoverGroup },
    });
  }, [configStore.hoverGroup]);

  useEffect(() => {
    setState({
      ...state,
      hoverLocation: { location: configStore.hoverLocation },
    });
  }, [configStore.hoverLocation]);

  useEffect(() => {
    setState({
      ...state,
      data: {
        ...state.data,
        selectedLocations: toJS(configStore.focusedLocations),
      },
    });
  }, [configStore.focusedLocations]);

  useEffect(() => {
    setState({
      ...state,
      data: {
        ...state.data,
        selectedGroups: toJS(configStore.selectedGroups),
      },
    });
  }, [configStore.selectedGroups]);

  useEffect(() => {
    if (UIStore.caseDataState !== ASYNC_STATES.SUCCEEDED) {
      return;
    }

    setState({
      ...state,
      data: {
        ...state.data,
        location_by_group: processLocationByGroup(),
        selectedGroups: toJS(configStore.selectedGroups),
      },
    });
  }, [UIStore.caseDataState, plotSettingsStore.locationGroupHideReference]);

  useEffect(() => {
    setState({
      ...state,
      data: {
        location_by_group: processLocationByGroup(),
        selectedGroups: toJS(configStore.selectedGroups),
        selectedLocations: [],
      },
    });
  }, []);

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
  if (Object.keys(config.group_cols).includes(configStore.groupKey)) {
    xLabel += `${config.group_cols[configStore.groupKey].title} `;
  } else if (configStore.groupKey === GROUP_SNV) {
    if (configStore.dnaOrAa === DNA_OR_AA.DNA) {
      xLabel += 'NT';
    } else {
      xLabel += 'AA';
    }
    xLabel += ' SNV ';
  }
  xLabel += ' (Cumulative, All Sequences)';

  if (configStore.groupKey === GROUP_SNV) {
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
      <PlotOptions>
        {configStore.groupKey === GROUP_SNV && (
          <OptionCheckboxContainer>
            <label>
              <input
                type="checkbox"
                checked={plotSettingsStore.locationGroupHideReference}
                onChange={onChangeHideReference}
              />
              Hide Reference Group
            </label>
          </OptionCheckboxContainer>
        )}
        <div className="spacer" />
      </PlotOptions>
      <div style={{ width: `${width}` }}>
        <VegaEmbed
          ref={vegaRef}
          data={state.data}
          spec={state.spec}
          signalListeners={state.signalListeners}
          dataListeners={state.dataListeners}
          signals={{
            hoverLocation: state.hoverLocation,
            hoverGroup: state.hoverGroup,
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

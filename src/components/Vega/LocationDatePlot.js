import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import _ from 'underscore';
import { aggregate } from '../../utils/transform';
import {
  NORM_MODES,
  COUNT_MODES,
  DATE_BINS,
} from '../../constants/plotSettings';

import EmptyPlot from '../Common/EmptyPlot';
import VegaEmbed from '../../react_vega/VegaEmbed';
import WarningBox from '../Common/WarningBox';
import initialSpec from '../../vega_specs/location_date.vg.json';

const PlotContainer = styled.div``;

const WarningContainer = styled.div`
  display: ${({ show }) => (show ? 'flex' : 'none')};
  flex-direction: column;
  align-items: stretch;
  justify-content: flex-start;

  // colors from Bootstrap
  background-color: #fff3cd;
  border: 1px solid #aaa;
  border-radius: 5px;

  margin: 0px 10px;
  margin-bottom: 15px;
  padding: 10px 20px;

  .warning-header {
    display: flex;
    flex-direction: row;
    align-items: center;
    margin-bottom: 5px;
    .warning-title {
      font-size: 1.25em;
    }
    .spacer {
      flex-grow: 1;
    }
  }

  .warning-text {
    margin: 0px;
    font-weight: normal;
  }
`;
WarningContainer.defaultProps = {
  show: true,
};

const PlotOptions = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-start;
  margin-bottom: 10px;
  padding-right: 10px;

  .plot-title {
    font-size: 1.25em;
    margin-right: 10px;
    padding-right: 10px;
    padding-left: 18px;

    border-right: 1px solid #ccc;
  }

  .spacer {
    flex-grow: 1;
  }
`;

const SelectContainer = styled.div`
  margin-right: 8px;
  font-weight: normal;
  select {
    margin-left: 0.65em;
    padding: 1px 4px;
    border-radius: 3px;
  }
`;

const LocationDatePlot = observer(({ width }) => {
  const vegaRef = useRef();
  const { dataStore, configStore, plotSettingsStore } = useStores();

  const handleHoverLocation = (...args) => {
    // Don't fire the action if there's no change
    let hoverLocation = args[1] === null ? null : args[1]['location'];
    if (hoverLocation === configStore.hoverLocation) {
      return;
    }
    configStore.updateHoverLocation(hoverLocation);
  };

  const handleSelected = (...args) => {
    // console.log(args);
    // Don't fire if the selection is the same
    if (_.isEqual(args[1], configStore.focusedLocations)) {
      return;
    }
    configStore.updateFocusedLocations(args[1]);
  };

  const [state, setState] = useState({
    showWarning: true,
    data: {
      location_data: [],
      selectedGroups: [],
      selected: JSON.parse(JSON.stringify(configStore.focusedLocations)),
    },
    spec: JSON.parse(JSON.stringify(initialSpec)),
    signalListeners: {
      hoverLocation: _.throttle(handleHoverLocation, 100),
    },
    dataListeners: {
      selected: handleSelected,
    },
  });

  const onDismissWarning = () => {
    setState({
      ...state,
      showWarning: false,
    });
  };

  const onChangeNormMode = (event) =>
    plotSettingsStore.setLocationDateNormMode(event.target.value);
  const onChangeCountMode = (event) =>
    plotSettingsStore.setLocationDateCountMode(event.target.value);
  const onChangeDateBin = (event) => {
    plotSettingsStore.setLocationDateDateBin(event.target.value);
  };

  // There's some weird data persistence things going on with the dateBin
  // in this plot... so instead of passing it as a signal, just modify the
  // spec and force a hard re-render
  useEffect(() => {
    // setState({ ...state, dateBin: event.target.value })
    const spec = JSON.parse(JSON.stringify(initialSpec));
    // Set the date bin
    let dateBin;
    if (plotSettingsStore.locationDateDateBin === DATE_BINS.DATE_BIN_DAY) {
      dateBin = 1000 * 60 * 60 * 24;
    } else if (
      plotSettingsStore.locationDateDateBin === DATE_BINS.DATE_BIN_WEEK
    ) {
      dateBin = 1000 * 60 * 60 * 24 * 7;
    } else if (
      plotSettingsStore.locationDateDateBin === DATE_BINS.DATE_BIN_MONTH
    ) {
      dateBin = 1000 * 60 * 60 * 24 * 30;
    }

    const dateBinSignal = _.findWhere(spec.signals, { name: 'dateBin' });
    dateBinSignal['value'] = dateBin;

    setState({
      ...state,
      spec,
    });
  }, [plotSettingsStore.locationDateDateBin]);

  useEffect(() => {
    setState({
      ...state,
      data: {
        ...state.data,
        selected: JSON.parse(JSON.stringify(configStore.focusedLocations)),
      },
    });
  }, [configStore.focusedLocations]);

  useEffect(() => {
    let locationData = JSON.parse(JSON.stringify(dataStore.aggLocationData));

    locationData.forEach((row) => {
      if (!dataStore.groupsToKeep.includes(row.group)) {
        row.group = 'other';
      }
    });

    locationData = aggregate({
      data: locationData,
      groupby: ['location', 'date', 'group'],
      fields: ['cases_sum'],
      ops: ['sum'],
      as: ['cases_sum'],
    });

    setState({
      ...state,
      data: {
        ...state.data,
        location_data: locationData,
        selectedGroups: JSON.parse(JSON.stringify(configStore.selectedGroups)),
      },
    });
  }, [
    dataStore.aggLocationData,
    configStore.selectedGroups,
    dataStore.groupsToKeep,
  ]);

  // Set the normalization mode
  const yField =
    plotSettingsStore.locationDateNormMode === NORM_MODES.NORM_PERCENTAGES
      ? 'cases_norm'
      : 'cases_sum_agg';
  // Set cumulative mode
  const cumulativeWindow =
    plotSettingsStore.locationDateCountMode === COUNT_MODES.COUNT_NEW
      ? [0, 0]
      : [null, 0];

  let yLabel = '';
  if (
    plotSettingsStore.locationDateCountMode === COUNT_MODES.COUNT_CUMULATIVE
  ) {
    yLabel += 'Cumulative ';
  }
  if (plotSettingsStore.locationDateNormMode === NORM_MODES.NORM_PERCENTAGES) {
    yLabel += '% ';
  }
  yLabel += 'Sequences by ' + configStore.getGroupLabel();

  let plotTitle = '';
  if (
    plotSettingsStore.locationDateCountMode === COUNT_MODES.COUNT_CUMULATIVE
  ) {
    plotTitle += 'Cumulative ';
  } else if (
    plotSettingsStore.locationDateCountMode === COUNT_MODES.COUNT_NEW
  ) {
    plotTitle += 'New ';
  }
  plotTitle += configStore.getGroupLabel();
  plotTitle +=
    plotSettingsStore.locationDateNormMode === NORM_MODES.NORM_PERCENTAGES
      ? ' Percentages'
      : ' Counts';

  if (plotSettingsStore.locationDateDateBin === DATE_BINS.DATE_BIN_DAY) {
    plotTitle += ' by Day';
  } else if (
    plotSettingsStore.locationDateDateBin === DATE_BINS.DATE_BIN_WEEK
  ) {
    plotTitle += ' by Week';
  } else if (
    plotSettingsStore.locationDateDateBin === DATE_BINS.DATE_BIN_MONTH
  ) {
    plotTitle += ' by Month';
  }
  if (configStore.selectedGroups.length > 0) {
    plotTitle += ' (Selected ' + configStore.getGroupLabel() + 's Only)';
  }

  const renderPlot = () => {
    if (configStore.selectedLocationIds.length == 0) {
      return (
        <EmptyPlot height={200}>
          <p>
            No locations selected. Please select one or more locations from the
            sidebar, under &quot;Selected Locations&quot;, to compare counts of{' '}
            <b>{configStore.getGroupLabel()}</b> between them.
          </p>
        </EmptyPlot>
      );
    }

    if (configStore.selectedGroups.length == 0) {
      return (
        <EmptyPlot height={200}>
          <p>
            Please select a <b>{configStore.getGroupLabel()}</b> from the legend
            (above) or the group plot (below), to compare its counts between the
            selected locations.
          </p>
        </EmptyPlot>
      );
    }

    return (
      <>
        <PlotOptions>
          <span className="plot-title">{plotTitle}</span>
          <SelectContainer>
            <label>
              <select
                value={plotSettingsStore.locationDateCountMode}
                onChange={onChangeCountMode}
              >
                <option value={COUNT_MODES.COUNT_NEW}>New</option>
                <option value={COUNT_MODES.COUNT_CUMULATIVE}>Cumulative</option>
              </select>
            </label>
          </SelectContainer>
          sequences, shown as{' '}
          <SelectContainer>
            <label>
              <select
                value={plotSettingsStore.locationDateNormMode}
                onChange={onChangeNormMode}
              >
                <option value={NORM_MODES.NORM_COUNTS}>Counts</option>
                <option value={NORM_MODES.NORM_PERCENTAGES}>Percentages</option>
              </select>
            </label>
          </SelectContainer>
          grouped by{' '}
          <SelectContainer>
            <label>
              <select
                value={plotSettingsStore.locationDateDateBin}
                onChange={onChangeDateBin}
              >
                <option value={DATE_BINS.DATE_BIN_DAY}>Day</option>
                <option value={DATE_BINS.DATE_BIN_WEEK}>Week</option>
                <option value={DATE_BINS.DATE_BIN_MONTH}>Month</option>
              </select>
            </label>
          </SelectContainer>
        </PlotOptions>
        <div style={{ width: `${width}` }}>
          <VegaEmbed
            ref={vegaRef}
            data={state.data}
            recreateOnDatasets={['selectedGroups']}
            spec={state.spec}
            signalListeners={state.signalListeners}
            dataListeners={state.dataListeners}
            width={width}
            signals={{
              hoverLocation: { location: configStore.hoverLocation },
              yField,
              cumulativeWindow,
              yLabel,
            }}
            actions={false}
          />
        </div>
      </>
    );
  };

  return (
    <PlotContainer>
      <WarningBox
        show={state.showWarning}
        onDismiss={onDismissWarning}
        text="Inconsistent sampling in the underlying data can result in missing
          data and artefacts in this visualization. Please interpret this data
          with care."
      />
      {renderPlot()}
    </PlotContainer>
  );
});
LocationDatePlot.propTypes = {
  width: PropTypes.number,
};
LocationDatePlot.defaultProps = {
  width: 100,
};

export default LocationDatePlot;

import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import _ from 'underscore';

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
  const { dataStore } = useStores();

  const handleHoverLocation = (...args) => {
    // Don't fire the action if there's no change
    let hoverLocation = args[1] === null ? null : args[1]['location'];
    if (hoverLocation === dataStore.hoverLocation) {
      return;
    }
    dataStore.updateHoverLocation(hoverLocation);
  };

  const handleSelected = (...args) => {
    // console.log(args);
    // Don't fire if the selection is the same
    if (_.isEqual(args[1], dataStore.focusedLocations)) {
      return;
    }
    dataStore.updateFocusedLocations(args[1]);
  };

  const [state, setState] = useState({
    showWarning: true,
    data: {
      location_data: [],
      selectedGroups: [],
      selected: JSON.parse(JSON.stringify(dataStore.focusedLocations)),
    },
    spec: JSON.parse(JSON.stringify(initialSpec)),
    signalListeners: {
      hoverLocation: _.throttle(handleHoverLocation, 100),
    },
    dataListeners: {
      selected: handleSelected,
    },
    normMode: 'counts', // 'percentages' or 'counts'
    countMode: 'new', // 'new' or 'cumulative'
    dateBin: 'day', // 'day', 'week', 'month'
  });

  const onDismissWarning = () => {
    setState({
      ...state,
      showWarning: false,
    });
  };

  const onChangeNormMode = (event) =>
    setState({ ...state, normMode: event.target.value });
  const onChangeCountMode = (event) =>
    setState({ ...state, countMode: event.target.value });

  // There's some weird data persistence things going on with the dateBin
  // in this plot... so instead of passing it as a signal, just modify the
  // spec and force a hard re-render
  const onChangeDateBin = (event) => {
    // setState({ ...state, dateBin: event.target.value })
    const spec = JSON.parse(JSON.stringify(initialSpec));
    // Set the date bin
    let dateBin;
    if (event.target.value === 'day') {
      dateBin = 1000 * 60 * 60 * 24;
    } else if (event.target.value === 'week') {
      dateBin = 1000 * 60 * 60 * 24 * 7;
    } else if (event.target.value === 'month') {
      dateBin = 1000 * 60 * 60 * 24 * 30;
    }

    const dateBinSignal = _.findWhere(spec.signals, { name: 'dateBin' });
    dateBinSignal['value'] = dateBin;

    setState({
      ...state,
      dateBin: event.target.value,
      spec,
    });
  };

  useEffect(() => {
    setState({
      ...state,
      data: {
        ...state.data,
        selected: JSON.parse(JSON.stringify(dataStore.focusedLocations)),
      },
    });
  }, [dataStore.focusedLocations]);

  useEffect(() => {
    setState({
      ...state,
      data: {
        ...state.data,
        location_data: JSON.parse(JSON.stringify(dataStore.aggLocationData)),
        selectedGroups: JSON.parse(JSON.stringify(dataStore.selectedGroups)),
      },
    });
  }, [dataStore.aggLocationData, dataStore.selectedGroups]);

  // Set the normalization mode
  const yField =
    state.normMode === 'percentages' ? 'cases_norm' : 'cases_sum_agg';
  // Set cumulative mode
  const cumulativeWindow = state.countMode === 'new' ? [0, 0] : [null, 0];

  let yLabel = '';
  if (state.countMode === 'cumulative') {
    yLabel += 'Cumulative ';
  }
  if (state.normMode === 'percentages') {
    yLabel += '% ';
  }
  yLabel += 'Sequences by ' + dataStore.getGroupLabel();

  let plotTitle = '';
  if (state.countMode === 'cumulative') {
    plotTitle += 'Cumulative ';
  } else if (state.countMode === 'new') {
    plotTitle += 'New ';
  }
  plotTitle += dataStore.getGroupLabel();
  plotTitle += state.normMode === 'percentages' ? ' Percentages' : ' Counts';

  if (state.dateBin === 'day') {
    plotTitle += ' by Day';
  } else if (state.dateBin === 'week') {
    plotTitle += ' by Week';
  } else if (state.dateBin === 'month') {
    plotTitle += ' by Month';
  }
  if (dataStore.selectedGroups.length > 0) {
    plotTitle += ' (Selected ' + dataStore.getGroupLabel() + 's Only)';
  }

  const renderPlot = () => {
    if (dataStore.selectedLocationIds.length == 0) {
      return (
        <EmptyPlot height={200}>
          <p>
            No locations selected. Please select one or more locations from the
            sidebar, under &quot;Selected Locations&quot;, to compare counts of{' '}
            <b>{dataStore.getGroupLabel()}</b> between them.
          </p>
        </EmptyPlot>
      );
    }

    if (dataStore.selectedGroups.length == 0) {
      return (
        <EmptyPlot height={200}>
          <p>
            Please select a <b>{dataStore.getGroupLabel()}</b> from the legend
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
              <select value={state.countMode} onChange={onChangeCountMode}>
                <option value="new">New</option>
                <option value="cumulative">Cumulative</option>
              </select>
            </label>
          </SelectContainer>
          sequences, shown as{' '}
          <SelectContainer>
            <label>
              <select value={state.normMode} onChange={onChangeNormMode}>
                <option value="counts">Counts</option>
                <option value="percentages">Percentages</option>
              </select>
            </label>
          </SelectContainer>
          grouped by{' '}
          <SelectContainer>
            <label>
              <select value={state.dateBin} onChange={onChangeDateBin}>
                <option value="day">Day</option>
                <option value="week">Week</option>
                <option value="month">Month</option>
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
              hoverLocation: { location: dataStore.hoverLocation },
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

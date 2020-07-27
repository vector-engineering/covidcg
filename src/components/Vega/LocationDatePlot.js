import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import _ from 'underscore';

import VegaEmbed from '../../react_vega/VegaEmbed';
import initialSpec from '../../vega/location_date.vg.json';

const PlotContainer = styled.div``;

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
  const { covidStore } = useStores();

  const handleHoverLocation = (...args) => {
    // Don't fire the action if there's no change
    let hoverLocation = args[1] === null ? null : args[1]['location'];
    if (hoverLocation === covidStore.hoverLocation) {
      return;
    }
    covidStore.updateHoverLocation(hoverLocation);
  };

  const handleSelected = (...args) => {
    // console.log(args);
    // Don't fire if the selection is the same
    if (_.isEqual(args[1], covidStore.focusedLocations)) {
      return;
    }
    covidStore.updateFocusedLocations(args[1]);
  };

  const [state, setState] = useState({
    data: {
      location_data: [],
      selectedGroups: [],
      selected: JSON.parse(JSON.stringify(covidStore.focusedLocations)),
    },
    spec: JSON.parse(JSON.stringify(initialSpec)),
    signalListeners: {
      hoverLocation: _.throttle(handleHoverLocation, 100),
    },
    dataListeners: {
      selected: handleSelected,
    },
    hoverLocation: {},
    normMode: 'counts', // 'percentages' or 'counts'
    countMode: 'new', // 'new' or 'cumulative'
    dateBin: 'day', // 'day', 'week', 'month'
  });

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
        selected: JSON.parse(JSON.stringify(covidStore.focusedLocations)),
      },
    });
  }, [covidStore.focusedLocations]);

  useEffect(() => {
    setState({
      ...state,
      data: {
        ...state.data,
        location_data: JSON.parse(JSON.stringify(covidStore.aggLocationData)),
        selectedGroups: JSON.parse(JSON.stringify(covidStore.selectedGroups)),
      },
    });
  }, [covidStore.aggLocationData, covidStore.selectedGroups]);

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
  if (covidStore.groupKey === 'lineage') {
    yLabel += 'Sequences by Lineage';
  } else if (covidStore.groupKey === 'clade') {
    yLabel += 'Sequences by Clade';
  } else if (covidStore.groupKey === 'snp') {
    yLabel +=
      'Sequences by ' + (covidStore.dnaOrAa === 'dna' ? 'NT' : 'AA') + ' SNP';
  }

  let plotTitle = '';
  if (state.countMode === 'cumulative') {
    plotTitle += 'Cumulative ';
  } else if (state.countMode === 'new') {
    plotTitle += 'New ';
  }

  if (covidStore.groupKey === 'lineage') {
    plotTitle += 'Lineage ';
  } else if (covidStore.groupKey === 'clade') {
    plotTitle += 'Clade ';
  } else if (covidStore.groupKey === 'snp') {
    if (covidStore.dnaOrAa === 'dna') {
      plotTitle += 'NT';
    } else {
      plotTitle += 'AA';
    }
    plotTitle += ' SNP ';
  }
  plotTitle += state.normMode === 'percentages' ? 'Percentages' : 'Counts';

  if (state.dateBin === 'day') {
    plotTitle += ' by Day';
  } else if (state.dateBin === 'week') {
    plotTitle += ' by Week';
  } else if (state.dateBin === 'month') {
    plotTitle += ' by Month';
  }

  return (
    <PlotContainer>
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
            hoverLocation: { location: covidStore.hoverLocation },
            yField,
            cumulativeWindow,
            yLabel,
          }}
          actions={false}
        />
      </div>
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

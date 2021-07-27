import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { toJS } from 'mobx';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { aggregate } from '../../utils/transform';
import { throttle } from '../../utils/func';

import {
  NORM_MODES,
  COUNT_MODES,
  DATE_BINS,
  PLOT_DOWNLOAD_OPTIONS,
  GROUPS,
  GROUP_SNV,
  ASYNC_STATES,
} from '../../constants/defs.json';

import EmptyPlot from '../Common/EmptyPlot';
import VegaEmbed from '../../react_vega/VegaEmbed';
import WarningBox from '../Common/WarningBox';
import DropdownButton from '../Buttons/DropdownButton';
import SkeletonElement from '../Common/SkeletonElement';
import LoadingSpinner from '../Common/LoadingSpinner';
import { PlotOptions, OptionSelectContainer } from './Plot.styles';

import initialSpec from '../../vega_specs/location_date.vg.json';
import { formatSnv } from '../../utils/snpUtils';

const PlotContainer = styled.div``;

const LocationDatePlot = observer(({ width }) => {
  const vegaRef = useRef();
  const { dataStore, configStore, UIStore, plotSettingsStore, groupDataStore } =
    useStores();

  const handleHoverLocation = (...args) => {
    // Don't fire the action if there's no change
    let hoverLocation = args[1] === null ? null : args[1]['location'];
    if (hoverLocation === configStore.hoverLocation) {
      return;
    }
    configStore.updateHoverLocation(hoverLocation);
  };

  const handleSelected = (...args) => {
    configStore.updateFocusedLocations(args[1]);
  };

  const processLocationData = () => {
    console.log('PROCESS LOCATION DATE DATA');
    let locationData;
    if (configStore.groupKey === GROUP_SNV) {
      if (dataStore.aggLocationSelectedSnvsDate === undefined) {
        return [];
      }

      locationData = toJS(dataStore.aggLocationSelectedSnvsDate);
    } else {
      if (dataStore.aggLocationGroupDate === undefined) {
        return [];
      }

      locationData = toJS(dataStore.aggLocationGroupDate).map((record) => {
        record.color = groupDataStore.getGroupColor(
          configStore.groupKey,
          record.group_id
        );
        record.group = record.group_id;
        record.group_name = record.group_id;
        return record;
      });
    }

    if (configStore.groupKey === GROUP_SNV) {
      // Filter out 'All Other Sequences' group
      locationData = locationData.filter((row) => {
        return row.group !== GROUPS.ALL_OTHER_GROUP;
      });
    }

    locationData = aggregate({
      data: locationData,
      groupby: ['location', 'collection_date', 'group'],
      fields: ['counts', 'group_name'],
      ops: ['sum', 'first'],
      as: ['counts', 'group_name'],
    }).map((record) => {
      // Add location counts
      record.location_counts = dataStore.countsPerLocationMap[record.location];
      record.location_date_count = dataStore.countsPerLocationDateMap
        .get(record.location)
        .get(record.collection_date);
      record.cumulative_location_date_count =
        dataStore.cumulativeCountsPerLocationDateMap
          .get(record.location)
          .get(record.collection_date);
      return record;
    });

    // console.log(locationData);
    // console.log(JSON.stringify(locationData));

    return locationData;
  };

  const processSelectedGroups = () => {
    return JSON.parse(JSON.stringify(configStore.selectedGroups));
  };

  const processFocusedLocations = () => {
    return JSON.parse(JSON.stringify(configStore.focusedLocations));
  };

  // There's some weird data persistence things going on with the dateBin
  // in this plot... so instead of passing it as a signal, just modify the
  // spec and force a hard re-render
  const injectDateBinIntoSpec = () => {
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

    const dateBinSignal = spec.signals.find(
      (signal) => signal.name === 'dateBin'
    );
    dateBinSignal['value'] = dateBin;

    return spec;
  };

  const [state, setState] = useState({
    showWarning: true,
    data: {
      location_data: [],
      selectedGroups: [],
      selected: [],
    },
    hoverLocation: null,
    spec: injectDateBinIntoSpec(),
    signalListeners: {
      hoverLocation: throttle(handleHoverLocation, 100),
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

  const handleDownloadSelect = (option) => {
    // console.log(option);
    // TODO: use the plot options and configStore options to build a more descriptive filename
    //       something like new_lineages_by_day_S_2020-05-03-2020-05-15_NYC.png...
    if (option === PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_DATA) {
      dataStore.downloadAggLocationGroupDate();
    } else if (option === PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_PNG) {
      vegaRef.current.downloadImage('png', 'vega-export.png', 1);
    } else if (option === PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_PNG_2X) {
      vegaRef.current.downloadImage('png', 'vega-export.png', 2);
    } else if (option === PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_PNG_4X) {
      vegaRef.current.downloadImage('png', 'vega-export.png', 4);
    } else if (option === PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_SVG) {
      vegaRef.current.downloadImage('svg', 'vega-export.svg');
    }
  };

  // Trigger date bin spec injection when the date bin in the store changes
  useEffect(() => {
    const spec = injectDateBinIntoSpec();
    setState({
      ...state,
      spec,
    });
  }, [plotSettingsStore.locationDateDateBin]);

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
        selected: processFocusedLocations(),
      },
    });
  }, [configStore.focusedLocations]);

  const refreshData = () => {
    if (
      configStore.groupKey !== GROUP_SNV ||
      UIStore.snvDataState !== ASYNC_STATES.SUCCEEDED
    ) {
      return;
    }

    setState({
      ...state,
      data: {
        ...state.data,
        location_data: processLocationData(),
        selectedGroups: processSelectedGroups(),
      },
    });
  };

  // Refresh data on mount (i.e., tab change) or when data state changes
  useEffect(refreshData, [
    UIStore.caseDataState,
    UIStore.snvDataState,
    configStore.selectedGroups,
  ]);
  useEffect(refreshData, []);

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
        <SkeletonElement delay={2} height={400}>
          <LoadingSpinner />
        </SkeletonElement>
      </div>
    );
  }

  if (configStore.selectedLocationNodes.length == 0) {
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

  if (
    configStore.groupKey === GROUP_SNV &&
    state.data.location_data.length === 0
  ) {
    return (
      <EmptyPlot height={200}>
        <p>
          No sequences with {configStore.getGroupLabel()}s:{' '}
          {`${configStore.selectedGroups
            .map((group) => group.group)
            .join(' & ')}`}
        </p>
      </EmptyPlot>
    );
  }

  let yLabel = '';
  if (
    plotSettingsStore.locationDateCountMode === COUNT_MODES.COUNT_CUMULATIVE
  ) {
    yLabel += 'Cumulative ';
  }
  if (plotSettingsStore.locationDateNormMode === NORM_MODES.NORM_PERCENTAGES) {
    yLabel += '% ';
  }
  if (configStore.groupKey === GROUP_SNV) {
    yLabel += 'Sequences with this ' + configStore.getGroupLabel();
  } else {
    yLabel += 'Sequences by ' + configStore.getGroupLabel();
  }

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
    if (configStore.groupKey === GROUP_SNV) {
      plotTitle += ` (${configStore.selectedGroups
        .map((group) => formatSnv(group.group, configStore.dnaOrAa))
        .join(' & ')})`;
    } else {
      plotTitle += ` (${configStore.selectedGroups
        .map((group) => group.group)
        .join(', ')})`;
    }
  }

  return (
    <PlotContainer>
      <WarningBox show={state.showWarning} onDismiss={onDismissWarning}>
        Inconsistent sampling in the underlying data can result in missing data
        and artefacts in this visualization. Please interpret this data with
        care.
      </WarningBox>
      <PlotOptions>
        <span className="plot-title">{plotTitle}</span>
        <OptionSelectContainer>
          <label>
            <select
              value={plotSettingsStore.locationDateCountMode}
              onChange={onChangeCountMode}
            >
              <option value={COUNT_MODES.COUNT_NEW}>New</option>
              <option value={COUNT_MODES.COUNT_CUMULATIVE}>Cumulative</option>
            </select>
          </label>
        </OptionSelectContainer>
        sequences, shown as{' '}
        <OptionSelectContainer>
          <label>
            <select
              value={plotSettingsStore.locationDateNormMode}
              onChange={onChangeNormMode}
            >
              <option value={NORM_MODES.NORM_COUNTS}>Counts</option>
              <option value={NORM_MODES.NORM_PERCENTAGES}>Percentages</option>
            </select>
          </label>
        </OptionSelectContainer>
        grouped by{' '}
        <OptionSelectContainer>
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
        </OptionSelectContainer>
        <div className="spacer"></div>
        <DropdownButton
          text={'Download'}
          options={[
            PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_DATA,
            PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_PNG,
            PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_PNG_2X,
            PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_PNG_4X,
            PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_SVG,
          ]}
          onSelect={handleDownloadSelect}
        />
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
            percentages:
              plotSettingsStore.locationDateNormMode ===
              NORM_MODES.NORM_PERCENTAGES,
            cumulative:
              plotSettingsStore.locationDateCountMode ===
              COUNT_MODES.COUNT_CUMULATIVE,
            skipFiltering: configStore.groupKey === GROUP_SNV,
            hoverLocation: state.hoverLocation,
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

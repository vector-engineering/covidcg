import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { aggregate } from '../../utils/transform';
import { formatSnv } from '../../utils/snpData';
import _ from 'underscore';

import {
  NORM_MODES,
  COUNT_MODES,
  DATE_BINS,
} from '../../constants/plotSettings';
import { PLOT_DOWNLOAD_OPTIONS } from '../../constants/download';
import { GROUPS } from '../../constants/groups';
import { GROUP_KEYS } from '../../constants/config';
import { ASYNC_STATES } from '../../constants/UI';

import EmptyPlot from '../Common/EmptyPlot';
import VegaEmbed from '../../react_vega/VegaEmbed';
import WarningBox from '../Common/WarningBox';
import DropdownButton from '../Buttons/DropdownButton';
import SkeletonElement from '../Common/SkeletonElement';
import LoadingSpinner from '../Common/LoadingSpinner';
import { PlotOptions, OptionSelectContainer } from './Plot.styles';

import initialSpec from '../../vega_specs/location_date.vg.json';

const PlotContainer = styled.div``;

const LocationDatePlot = observer(({ width }) => {
  const vegaRef = useRef();
  const { dataStore, configStore, UIStore, plotSettingsStore } = useStores();

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

  const processLocationData = () => {
    let locationData;
    if (configStore.groupKey === GROUP_KEYS.GROUP_SNV) {
      if (dataStore.dataAggLocationSnvDate === undefined) {
        return [];
      }

      locationData = JSON.parse(
        JSON.stringify(dataStore.dataAggLocationSnvDate)
      );
      // Filter out 'Other' group
      locationData = locationData.filter((row) => {
        return row.group !== GROUPS.ALL_OTHER_GROUP;
      });
    } else {
      if (dataStore.dataAggLocationGroupDate === undefined) {
        return [];
      }

      locationData = JSON.parse(
        JSON.stringify(dataStore.dataAggLocationGroupDate)
      );

      locationData.forEach((row) => {
        if (!dataStore.groupsToKeep.includes(row.group)) {
          row.group = GROUPS.OTHER_GROUP;
          row.groupName = GROUPS.OTHER_GROUP;
        }
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
      fields: ['cases_sum', 'location_counts'],
      ops: ['sum', 'max'],
      as: ['cases_sum', 'location_counts'],
    });

    // Manually join the countsPerLocationDate to locationData
    locationData.forEach((row) => {
      row.location_date_count =
        dataStore.countsPerLocationDate[row.location][row.date.toString()];
    });

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

    const dateBinSignal = _.findWhere(spec.signals, { name: 'dateBin' });
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
    spec: injectDateBinIntoSpec(),
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

  const handleDownloadSelect = (option) => {
    // console.log(option);
    // TODO: use the plot options and configStore options to build a more descriptive filename
    //       something like new_lineages_by_day_S_2020-05-03-2020-05-15_NYC.png...
    if (option === PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_DATA) {
      dataStore.downloadDataAggLocationGroupDate();
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
      data: {
        ...state.data,
        selected: processFocusedLocations(),
      },
    });
  }, [configStore.focusedLocations]);

  useEffect(() => {
    if (
      configStore.groupKey !== GROUP_KEYS.GROUP_SNV ||
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
  }, [
    UIStore.snvDataState,
    configStore.selectedGroups,
    configStore.dateRange,
    dataStore.groupsToKeep,
  ]);

  useEffect(() => {
    if (
      configStore.groupKey === GROUP_KEYS.GROUP_SNV ||
      UIStore.caseDataState !== ASYNC_STATES.SUCCEEDED
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
  }, [
    UIStore.caseDataState,
    configStore.selectedGroups,
    configStore.dateRange,
    dataStore.groupsToKeep,
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
    configStore.groupKey === GROUP_KEYS.GROUP_SNV &&
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
  if (configStore.groupKey === GROUP_KEYS.GROUP_SNV) {
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
    if (configStore.groupKey === GROUP_KEYS.GROUP_SNV) {
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
      <WarningBox
        show={state.showWarning}
        onDismiss={onDismissWarning}
        text="Inconsistent sampling in the underlying data can result in missing
          data and artefacts in this visualization. Please interpret this data
          with care."
      />
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
            skipFiltering: configStore.groupKey === GROUP_KEYS.GROUP_SNV,
            hoverLocation: { location: configStore.hoverLocation },
            yField,
            cumulativeWindow,
            yLabel,
            yFormat:
              plotSettingsStore.locationDateNormMode === NORM_MODES.NORM_COUNTS
                ? 's'
                : '%',
            tooltipCountFormat:
              plotSettingsStore.locationDateNormMode === NORM_MODES.NORM_COUNTS
                ? 'd'
                : '.1%',
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

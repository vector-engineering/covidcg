import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import { toJS } from 'mobx';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import {
  ASYNC_STATES,
  NORM_MODES,
  COUNT_MODES,
  DATE_BINS,
  GROUP_SNV,
  PLOT_DOWNLOAD_OPTIONS,
  GROUPS,
} from '../../constants/defs.json';
import { throttle } from '../../utils/func';

import EmptyPlot from '../Common/EmptyPlot';
import WarningBox from '../Common/WarningBox';
import DropdownButton from '../Buttons/DropdownButton';
import VegaEmbed from '../../react_vega/VegaEmbed';
import SkeletonElement from '../Common/SkeletonElement';
import { PlotTitle, PlotOptions, OptionSelectContainer } from './Plot.styles';

import initialSpec from '../../vega_specs/group_stack.vg.json';

const GroupStackPlot = observer(({ width }) => {
  const vegaRef = useRef();
  const { dataStore, UIStore, configStore, plotSettingsStore, groupDataStore } =
    useStores();

  // disable this for now
  // const handleBrush = (...args) => {
  //   let dateRange = args[1];

  //   if (dateRange === null) {
  //     // Reset time range
  //     configStore.updateDateRange([-1, -1]);
  //   } else if (
  //     dateRange[0] === configStore.dateRange[0] &&
  //     dateRange[1] === configStore.dateRange[1]
  //   ) {
  //     // No change, return
  //   } else if (dateRange !== null) {
  //     configStore.updateDateRange([
  //       dateRange[0].getTime(),
  //       dateRange[1].getTime(),
  //     ]);
  //   }
  // };

  const handleHoverGroup = (...args) => {
    // Don't fire the action if there's no change
    let hoverGroup = args[1] === null ? null : args[1]['group'];
    if (hoverGroup === configStore.hoverGroup) {
      return;
    }
    // Ignore for some special groups
    if (hoverGroup === GROUPS.ALL_OTHER_GROUP) {
      return;
    }
    configStore.updateHoverGroup(hoverGroup);
  };

  const handleSelected = (...args) => {
    // console.log(args);

    // Ignore selections in SNV mode
    if (configStore.groupKey === GROUP_SNV) {
      return;
    }

    configStore.updateSelectedGroups(args[1]);
  };

  const handleDownloadSelect = (option) => {
    // console.log(option);
    // TODO: use the plot options and configStore options to build a more descriptive filename
    //       something like new_lineages_by_day_S_2020-05-03-2020-05-15_NYC.png...
    if (option === PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_DATA) {
      dataStore.downloadAggGroupDate();
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

  const processData = () => {
    console.log('GROUP STACK PROCESS DATA');
    // console.log(dataStore.aggGroupDate);

    if (configStore.groupKey === GROUP_SNV) {
      return toJS(dataStore.aggSelectedSnvsDate);
    }

    // For non-SNV mode, we'll need some additional fields:
    // 1) color of group
    // 2) name of group (same as group id)
    return toJS(dataStore.aggGroupDate).map((record) => {
      record.color = groupDataStore.getGroupColor(
        configStore.groupKey,
        record.group_id
      );
      record.group = record.group_id;
      record.group_name = record.group_id;
      return record;
    });
  };

  const [state, setState] = useState({
    showWarning: true,
    // data: {},
    hoverGroup: null,
    signalListeners: {
      // detailDomain: debounce(handleBrush, 500),
      hoverBar: throttle(handleHoverGroup, 100),
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
    plotSettingsStore.setGroupStackNormMode(event.target.value);
  const onChangeCountMode = (event) =>
    plotSettingsStore.setGroupStackCountMode(event.target.value);
  const onChangeDateBin = (event) =>
    plotSettingsStore.setGroupStackDateBin(event.target.value);

  useEffect(() => {
    setState({
      ...state,
      hoverGroup: { group: configStore.hoverGroup },
    });
  }, [configStore.hoverGroup]);

  // Update internal selected groups copy
  useEffect(() => {
    // console.log('SELECTED GROUPS');
    // Skip this update if we're in SNV mode
    if (configStore.groupKey === GROUP_SNV) {
      return;
    }

    setState({
      ...state,
      data: {
        ...state.data,
        selected: toJS(configStore.selectedGroups),
      },
    });
  }, [configStore.selectedGroups]);

  const refreshData = () => {
    // Skip this if we're not in SNV mode
    if (configStore.groupKey !== GROUP_SNV) {
      return;
    }

    // Skip unless the SNV data finished processing
    if (UIStore.snvDataState !== ASYNC_STATES.SUCCEEDED) {
      return;
    }

    if (UIStore.caseDataState !== ASYNC_STATES.SUCCEEDED) {
      return;
    }

    setState({
      ...state,
      data: {
        cases_by_date_and_group: processData(),
        selected: toJS(configStore.selectedGroups),
      },
    });
  };

  // Refresh data on mount (i.e., tab change) or when data state changes
  useEffect(refreshData, [UIStore.caseDataState, UIStore.snvDataState]);
  useEffect(refreshData, []);

  // For development in Vega Editor
  // console.log(JSON.stringify(caseData));

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
        <SkeletonElement delay={2} height={400} />
      </div>
    );
  }

  if (configStore.selectedLocationNodes.length === 0) {
    return (
      <EmptyPlot height={250}>
        <p>
          No locations selected. Please select one or more locations from the
          sidebar, under &quot;Selected Locations&quot;, to compare counts of{' '}
          <b>{configStore.getGroupLabel()}</b> between them.
        </p>
      </EmptyPlot>
    );
  }

  let plotTitle = '';
  if (plotSettingsStore.groupStackCountMode === COUNT_MODES.COUNT_CUMULATIVE) {
    plotTitle += 'Cumulative ';
  } else if (plotSettingsStore.groupStackCountMode === COUNT_MODES.COUNT_NEW) {
    plotTitle += 'New ';
  }
  plotTitle += configStore.getGroupLabel();
  plotTitle +=
    plotSettingsStore.groupStackNormMode === NORM_MODES.NORM_PERCENTAGES
      ? ' Percentages'
      : ' Counts';

  if (plotSettingsStore.groupStackDateBin === DATE_BINS.DATE_BIN_DAY) {
    plotTitle += ' by Day';
  } else if (plotSettingsStore.groupStackDateBin === DATE_BINS.DATE_BIN_WEEK) {
    plotTitle += ' by Week';
  } else if (plotSettingsStore.groupStackDateBin === DATE_BINS.DATE_BIN_MONTH) {
    plotTitle += ' by Month';
  }

  const maxShownLocations = 4;
  let selectedLocationsText =
    '(' +
    configStore.selectedLocationNodes
      .slice(0, maxShownLocations)
      .map((node) => node.value)
      .join(', ');
  if (configStore.selectedLocationNodes.length > maxShownLocations) {
    selectedLocationsText += ', ...)';
  } else {
    selectedLocationsText += ')';
  }

  // Set the stack offset mode
  const stackOffset =
    plotSettingsStore.groupStackNormMode === NORM_MODES.NORM_PERCENTAGES
      ? 'normalize'
      : 'zero';
  // Set the date bin
  let dateBin;
  if (plotSettingsStore.groupStackDateBin === DATE_BINS.DATE_BIN_DAY) {
    dateBin = 1000 * 60 * 60 * 24;
  } else if (plotSettingsStore.groupStackDateBin === DATE_BINS.DATE_BIN_WEEK) {
    dateBin = 1000 * 60 * 60 * 24 * 7;
  } else if (plotSettingsStore.groupStackDateBin === DATE_BINS.DATE_BIN_MONTH) {
    dateBin = 1000 * 60 * 60 * 24 * 30;
  }
  // If running in cumulative mode, add the vega transformation
  // By default the cumulative transformation is dumped into a column
  // "counts_cumulative", so if active, just overwrite the "counts"
  // column with this cumulative count
  const cumulativeWindow =
    plotSettingsStore.groupStackCountMode === COUNT_MODES.COUNT_CUMULATIVE
      ? [null, 0]
      : [0, 0];

  // Adapt labels to groupings
  let detailYLabel = '';
  if (plotSettingsStore.groupStackCountMode === COUNT_MODES.COUNT_CUMULATIVE) {
    detailYLabel += 'Cumulative ';
  }
  if (plotSettingsStore.groupStackNormMode === NORM_MODES.NORM_PERCENTAGES) {
    detailYLabel += '% ';
  }
  detailYLabel += 'Sequences by ' + configStore.getGroupLabel();

  // Hide the detail view in SNV mode when there's no selections
  // Also disable the plot options when the detail panel is hidden
  const hideDetail =
    configStore.groupKey === GROUP_SNV &&
    configStore.selectedGroups.length === 0;
  const detailHeight = hideDetail ? 0 : 280;
  const height = hideDetail ? 60 : 380;

  return (
    <div>
      <WarningBox show={state.showWarning} onDismiss={onDismissWarning}>
        Inconsistent sampling in the underlying data can result in missing data
        and artefacts in this visualization. Please interpret this data with
        care.
      </WarningBox>
      {hideDetail && (
        <EmptyPlot height={100}>
          <p>
            No {configStore.getGroupLabel()}s selected. Please select one or
            more {configStore.getGroupLabel()}s from the legend, frequency plot,
            or table.
          </p>
        </EmptyPlot>
      )}
      {!hideDetail && (
        <PlotOptions>
          <PlotTitle>
            <span className="title">{plotTitle}</span>
            <span className="subtitle">{selectedLocationsText}</span>
          </PlotTitle>
          <OptionSelectContainer>
            <label>
              <select
                value={plotSettingsStore.groupStackCountMode}
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
                value={plotSettingsStore.groupStackNormMode}
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
                value={plotSettingsStore.groupStackDateBin}
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
      )}

      <div style={{ width: `${width}px` }}>
        <VegaEmbed
          ref={vegaRef}
          data={state.data}
          spec={initialSpec}
          signalListeners={state.signalListeners}
          dataListeners={state.dataListeners}
          signals={{
            disableSelectionColoring: configStore.groupKey === GROUP_SNV,
            detailHeight,
            hoverBar: state.hoverGroup,
            stackOffset,
            dateBin,
            cumulativeWindow,
            detailYLabel,
            yFormat:
              plotSettingsStore.groupStackNormMode === NORM_MODES.NORM_COUNTS
                ? 's'
                : '%',
          }}
          cheapSignals={['hoverBar']}
          width={width}
          height={height}
          actions={false}
        />
      </div>
    </div>
  );
});

GroupStackPlot.propTypes = {
  width: PropTypes.number.isRequired,
};

export default GroupStackPlot;

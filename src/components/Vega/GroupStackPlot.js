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
  GROUP_MUTATION,
  PLOT_DOWNLOAD_OPTIONS,
  GROUPS,
} from '../../constants/defs.json';
import { throttle } from '../../utils/func';
import { getValidGroups } from '../../utils/data';
import { aggregate } from '../../utils/transform';

import LowFreqFilter from './LowFreqFilter';
import EmptyPlot from '../Common/EmptyPlot';
import WarningBox from '../Common/WarningBox';
import DropdownButton from '../Buttons/DropdownButton';
import VegaEmbed from '../../react_vega/VegaEmbed';
import SkeletonElement from '../Common/SkeletonElement';
import StyledDropdownTreeSelect from '../Common/StyledDropdownTreeSelect';
import { PlotTitle, OptionSelectContainer } from './Plot.styles';
import { PlotHeader, PlotOptionsRow } from './GroupStackPlot.styles';

import initialSpec from '../../vega_specs/group_stack.vg.json';

const GroupStackPlot = observer(({ width }) => {
  const vegaRef = useRef();
  const { dataStore, UIStore, configStore, plotSettingsStore, groupDataStore } =
    useStores();

  const handleHoverGroup = (...args) => {
    // Don't fire the action if there's no change
    configStore.updateHoverGroup(args[1] === null ? null : args[1]['group']);
  };

  const handleSelected = (...args) => {
    // Ignore selections in mutation mode
    if (configStore.groupKey === GROUP_MUTATION) {
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
    // console.log('GROUP STACK PROCESS DATA');

    if (configStore.groupKey === GROUP_MUTATION) {
      let data = toJS(dataStore.aggLocationSelectedMutationsDate);

      // Filter focused locations
      const focusedLocations = state.focusLocationTree
        .filter((node) => node.checked)
        .map((node) => node.value);
      data = data.filter((record) =>
        focusedLocations.includes(record.location)
      );

      // Re-aggregate
      data = aggregate({
        data,
        groupby: ['group', 'collection_date'],
        fields: ['counts', 'color', 'group_name'],
        ops: ['sum', 'first', 'first'],
        as: ['counts', 'color', 'group_name'],
      });

      return data;
    }

    // For non-mutation mode, we'll need some additional fields:
    // 1) color of group
    // 2) name of group (same as group id)
    // Also collapse low-frequency groups based on settings
    const validGroups = getValidGroups({
      records: dataStore.groupCounts,
      lowFreqFilterType: plotSettingsStore.groupStackLowFreqFilter,
      lowFreqFilterValue: plotSettingsStore.groupStackLowFreqValue,
    });
    let data = toJS(dataStore.aggLocationGroupDate);
    
    // Filter focused locations
    const focusedLocations = state.focusLocationTree
      .filter((node) => node.checked)
      .map((node) => node.value);
    data = data.filter((record) =>
      focusedLocations.includes(record.location)
    );
    
    data = data.map((record) => {
      if (!validGroups.includes(record.group_id)) {
        record.group = GROUPS.OTHER_GROUP;
      } else {
        record.group = record.group_id;
      }

      record.group_name = record.group;
      record.color = groupDataStore.getGroupColor(
        configStore.groupKey,
        record.group
      );

      return record;
    });

    data = aggregate({
      data,
      groupby: ['group', 'collection_date'],
      fields: ['counts', 'color', 'group_name'],
      ops: ['sum', 'first', 'first'],
      as: ['counts', 'color', 'group_name'],
    });

    // console.log(data);

    return data;
  };

  const [state, setState] = useState({
    showWarning: true,
    // data: {},
    hoverGroup: null,
    focusLocationTree: [],
    signalListeners: {
      // detailDomain: debounce(handleBrush, 500),
      hoverBar: throttle(handleHoverGroup, 100),
    },
    dataListeners: {
      selected: handleSelected,
    },
  });

  // Update state based on the focused location dropdown select
  const focusLocationSelectOnChange = (event) => {
    let focusLocationTree = state.focusLocationTree.slice();
    focusLocationTree.forEach((node) => {
      if (node.value === event.value) {
        node.checked = event.checked;
      }
    });

    setState({
      ...state,
      focusLocationTree,
    });
  };

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
    // Skip this update if we're in mutation mode
    if (configStore.groupKey === GROUP_MUTATION) {
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
    // Skip unless the mutation data finished processing
    if (
      configStore.groupKey === GROUP_MUTATION &&
      UIStore.mutationDataState !== ASYNC_STATES.SUCCEEDED
    ) {
      return;
    }

    if (UIStore.caseDataState !== ASYNC_STATES.SUCCEEDED) {
      return;
    }

    // console.log('REFRESH DATA FROM STORE');

    // Update focus location tree with new locations
    const focusLocationTree = [];
    Object.keys(dataStore.countsPerLocationMap).forEach((loc) => {
      focusLocationTree.push({
        label: loc + ' (' + dataStore.countsPerLocationMap[loc] + ')',
        value: loc,
        checked: true,
      });
    });

    const newState = {
      ...state,
      data: {
        cases_by_date_and_group: processData(),
        selected: toJS(configStore.selectedGroups),
      },
    };

    // Only update location tree if the location list changed
    // This is to prevent firing the processData event twice when we change selected groups
    const prevLocs = state.focusLocationTree.map((node) => node.value);
    if (
      focusLocationTree.length !== state.focusLocationTree.length ||
      !focusLocationTree.every((node) => prevLocs.includes(node.value))
    ) {
      newState['focusLocationTree'] = focusLocationTree;
    }

    setState(newState);
  };

  useEffect(() => {
    // console.log('UPDATE DATA FROM FOCUSED LOCATIONS');

    setState({
      ...state,
      data: {
        ...state.data,
        cases_by_date_and_group: processData(),
      },
    });
  }, [state.focusLocationTree]);

  const focusLocationDropdownContainer = (
    <StyledDropdownTreeSelect
      mode={'multiSelect'}
      data={state.focusLocationTree}
      className="geo-dropdown-tree-select"
      clearSearchOnChange={false}
      keepTreeOnSearch={true}
      keepChildrenOnSearch={true}
      showPartiallySelected={false}
      inlineSearchInput={false}
      texts={{
        placeholder: 'Search...',
        noMatches: 'No matches found',
      }}
      onChange={focusLocationSelectOnChange}
      // onAction={treeSelectOnAction}
      // onNodeToggle={treeSelectOnNodeToggleCurrentNode}
    />
  );

  // Refresh data on mount (i.e., tab change) or when data state changes
  useEffect(refreshData, [
    UIStore.caseDataState,
    UIStore.mutationDataState,
    plotSettingsStore.groupStackLowFreqFilter,
    plotSettingsStore.groupStackLowFreqValue,
  ]);
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
    state.focusLocationTree
      .filter((node) => node.checked)
      .slice(0, maxShownLocations)
      .map((node) => node.value)
      .join(', ');
  if (state.focusLocationTree.length > maxShownLocations) {
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

  // Hide the detail view in mutation mode when there's no selections
  // Also disable the plot options when the detail panel is hidden
  const hideDetail =
    configStore.groupKey === GROUP_MUTATION &&
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
        <PlotHeader>
          <PlotTitle style={{ gridRow: '1/-1' }}>
            <span className="title">{plotTitle}</span>
            <span className="subtitle">{selectedLocationsText}</span>
          </PlotTitle>
          <PlotOptionsRow>
            <OptionSelectContainer>
              <label>
                <select
                  value={plotSettingsStore.groupStackCountMode}
                  onChange={onChangeCountMode}
                >
                  <option value={COUNT_MODES.COUNT_NEW}>New</option>
                  <option value={COUNT_MODES.COUNT_CUMULATIVE}>
                    Cumulative
                  </option>
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
                  <option value={NORM_MODES.NORM_PERCENTAGES}>
                    Percentages
                  </option>
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
          </PlotOptionsRow>
          {configStore.groupKey !== GROUP_MUTATION && (
            <PlotOptionsRow>
              <LowFreqFilter
                lowFreqFilterType={plotSettingsStore.groupStackLowFreqFilter}
                lowFreqFilterValue={plotSettingsStore.groupStackLowFreqValue}
                updateLowFreqFilterType={
                  plotSettingsStore.setGroupStackLowFreqFilter
                }
                updateLowFreqFilterValue={
                  plotSettingsStore.setGroupStackLowFreqValue
                }
              ></LowFreqFilter>
            </PlotOptionsRow>
          )}
          <PlotOptionsRow>Only show locations:&nbsp;{focusLocationDropdownContainer}</PlotOptionsRow>
          <PlotOptionsRow style={{ justifyContent: 'flex-end' }}>
            <DropdownButton
              text={'Download'}
              options={[
                PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_DATA,
                PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_PNG,
                PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_PNG_2X,
                PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_PNG_4X,
                PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_SVG,
              ]}
              style={{ minWidth: '90px' }}
              onSelect={handleDownloadSelect}
            />
          </PlotOptionsRow>
        </PlotHeader>
      )}

      <div style={{ width: `${width}px` }}>
        <VegaEmbed
          ref={vegaRef}
          data={state.data}
          spec={initialSpec}
          signalListeners={state.signalListeners}
          dataListeners={state.dataListeners}
          signals={{
            disableSelectionColoring: configStore.groupKey === GROUP_MUTATION,
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

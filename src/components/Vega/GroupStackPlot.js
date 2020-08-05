import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { toJS } from 'mobx';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { ASYNC_STATES } from '../../constants/UI';
import {
  NORM_MODES,
  COUNT_MODES,
  DATE_BINS,
} from '../../constants/plotSettings';
import { aggregate } from '../../utils/transform';
import _ from 'underscore';

import EmptyPlot from '../Common/EmptyPlot';
import WarningBox from '../Common/WarningBox';
import DropdownButton from '../Buttons/DropdownButton';
import VegaEmbed from '../../react_vega/VegaEmbed';
import SkeletonElement from '../Common/SkeletonElement';
import LoadingSpinner from '../Common/LoadingSpinner';

import { mergeGroupsIntoOther } from './utils';
import initialSpec from '../../vega_specs/bar_stack_v1.vg.json';

const PlotOptions = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-start;
  margin-bottom: 10px;
  padding-right: 10px;

  .area-stack-title {
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

const GroupStackPlot = observer(({ width }) => {
  const vegaRef = useRef();
  const { dataStore, UIStore, configStore, plotSettingsStore } = useStores();

  const handleBrush = (...args) => {
    let dateRange = args[1];
    if (dateRange !== null) {
      configStore.selectDateRange([
        dateRange[0].getTime(),
        dateRange[1].getTime(),
      ]);
    } else {
      // Reset time range
      configStore.selectDateRange([-1, -1]);
    }
  };

  const handleHoverGroup = (...args) => {
    // Don't fire the action if there's no change
    let hoverGroup = args[1] === null ? null : args[1]['group'];
    if (hoverGroup === configStore.hoverGroup) {
      return;
    }
    // console.log('Updating store hoverGroup from plot', hoverGroup);
    // console.log('state hovergroup', state.hoverGroup);
    configStore.updateHoverGroup(hoverGroup);
  };

  const handleSelected = (...args) => {
    // console.log(args);

    // Don't fire if the selection is the same
    if (_.isEqual(args[1], configStore.selectedGroups)) {
      return;
    }

    configStore.updateSelectedGroups(args[1]);
  };

  const handleDownloadSelect = (option) => {
    // console.log(option);
    // TODO: use the plot options and configStore options to build a more descriptive filename
    //       something like new_lineages_by_day_S_2020-05-03-2020-05-15_NYC.png...
    if (option === 'PNG') {
      vegaRef.current.downloadImage('png', 'vega-export.png', 1);
    } else if (option === 'PNG (2X)') {
      vegaRef.current.downloadImage('png', 'vega-export.png', 2);
    } else if (option === 'PNG (4X)') {
      vegaRef.current.downloadImage('png', 'vega-export.png', 4);
    } else if (option === 'SVG') {
      vegaRef.current.downloadImage('svg', 'vega-export.svg');
    }
  };

  const initialData = aggregate({
    data: toJS(dataStore.caseData),
    groupby: ['date', 'group'],
    fields: ['cases_sum', 'color'],
    ops: ['sum', 'max'],
    as: ['cases_sum', 'color'],
  });

  const [state, setState] = useState({
    showWarning: true,
    data: {
      cases_by_date_and_group: mergeGroupsIntoOther(
        initialData,
        dataStore.groupsToKeep
      ),
      selected: JSON.parse(JSON.stringify(configStore.selectedGroups)),
    },
    signalListeners: {
      detailDomain: _.debounce(handleBrush, 500),
      hoverBar: _.throttle(handleHoverGroup, 100),
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

  // Update internal caseData copy
  useEffect(() => {
    setState({
      ...state,
      data: {
        ...state.data,
        cases_by_date_and_group: mergeGroupsIntoOther(
          aggregate({
            data: toJS(dataStore.caseData),
            groupby: ['date', 'group'],
            fields: ['cases_sum', 'color'],
            ops: ['sum', 'max'],
            as: ['cases_sum', 'color'],
          }),
          dataStore.groupsToKeep
        ),
      },
    });
  }, [dataStore.caseData, dataStore.groupsToKeep]);

  // Update internal selected groups copy
  useEffect(() => {
    // console.log('Update selectedGroups');
    setState({
      ...state,
      data: {
        ...state.data,
        selected: JSON.parse(JSON.stringify(configStore.selectedGroups)),
      },
    });
  }, [configStore.selectedGroups]);

  // For development in Vega Editor
  // console.log(JSON.stringify(caseData));

  let areaStackTitle = '';
  if (plotSettingsStore.groupStackCountMode === COUNT_MODES.COUNT_CUMULATIVE) {
    areaStackTitle += 'Cumulative ';
  } else if (plotSettingsStore.groupStackCountMode === COUNT_MODES.COUNT_NEW) {
    areaStackTitle += 'New ';
  }
  areaStackTitle += configStore.getGroupLabel();
  areaStackTitle +=
    plotSettingsStore.groupStackNormMode === NORM_MODES.NORM_PERCENTAGES
      ? ' Percentages'
      : ' Counts';

  if (plotSettingsStore.groupStackDateBin === DATE_BINS.DATE_BIN_DAY) {
    areaStackTitle += ' by Day';
  } else if (plotSettingsStore.groupStackDateBin === DATE_BINS.DATE_BIN_WEEK) {
    areaStackTitle += ' by Week';
  } else if (plotSettingsStore.groupStackDateBin === DATE_BINS.DATE_BIN_MONTH) {
    areaStackTitle += ' by Month';
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
  // "cases_sum_cumulative", so if active, just overwrite the "cases_sum"
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

  if (configStore.selectedLocationIds.length === 0) {
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

  return (
    <div>
      <WarningBox
        show={state.showWarning}
        onDismiss={onDismissWarning}
        text="Inconsistent sampling in the underlying data can result in missing
          data and artefacts in this visualization. Please interpret this data
          with care."
      />
      <PlotOptions>
        <span className="area-stack-title">{areaStackTitle}</span>
        <SelectContainer>
          <label>
            <select
              value={plotSettingsStore.groupStackCountMode}
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
              value={plotSettingsStore.groupStackNormMode}
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
              value={plotSettingsStore.groupStackDateBin}
              onChange={onChangeDateBin}
            >
              <option value={DATE_BINS.DATE_BIN_DAY}>Day</option>
              <option value={DATE_BINS.DATE_BIN_WEEK}>Week</option>
              <option value={DATE_BINS.DATE_BIN_MONTH}>Month</option>
            </select>
          </label>
        </SelectContainer>
        <div className="spacer"></div>
        <DropdownButton
          text={'Download'}
          options={['PNG', 'PNG (2X)', 'PNG (4X)', 'SVG']}
          onSelect={handleDownloadSelect}
        />
      </PlotOptions>

      <div style={{ width: `${width}px` }}>
        <VegaEmbed
          ref={vegaRef}
          data={state.data}
          spec={initialSpec}
          signalListeners={state.signalListeners}
          dataListeners={state.dataListeners}
          signals={{
            hoverBar: { group: configStore.hoverGroup },
            stackOffset,
            dateBin,
            cumulativeWindow,
            detailYLabel,
          }}
          cheapSignals={['hoverBar']}
          width={width}
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

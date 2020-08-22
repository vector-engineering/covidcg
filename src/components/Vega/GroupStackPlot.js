import React, { useState, useEffect, useRef } from 'react';
import { toJS } from 'mobx';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { ASYNC_STATES } from '../../constants/UI';
import {
  NORM_MODES,
  COUNT_MODES,
  DATE_BINS,
} from '../../constants/plotSettings';
import { GROUP_KEYS } from '../../constants/config';
import { PLOT_DOWNLOAD_OPTIONS } from '../../constants/download';
import {
  dnaSnpToInt,
  geneAaSnpToInt,
  proteinAaSnpToInt,
} from '../../utils/snpData';
import _ from 'underscore';

import EmptyPlot from '../Common/EmptyPlot';
import WarningBox from '../Common/WarningBox';
import DropdownButton from '../Buttons/DropdownButton';
import VegaEmbed from '../../react_vega/VegaEmbed';
import SkeletonElement from '../Common/SkeletonElement';
import LoadingSpinner from '../Common/LoadingSpinner';
import { PlotTitle, PlotOptions, OptionSelectContainer } from './Plot.styles';

import { mergeGroupsIntoOther } from './utils';
import initialSpec from '../../vega_specs/bar_stack_v1.vg.json';

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
    if (option === PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_DATA) {
      dataStore.downloadDataAggGroupDate();
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

  const [state, setState] = useState({
    showWarning: true,
    data: {
      cases_by_date_and_group: mergeGroupsIntoOther(
        JSON.parse(JSON.stringify(dataStore.dataAggGroupDate)),
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
          JSON.parse(JSON.stringify(dataStore.dataAggGroupDate)),
          dataStore.groupsToKeep
        ),
      },
    });
  }, [dataStore.selectedRowsHash, dataStore.groupsToKeep]);

  // Update internal selected groups copy
  useEffect(() => {
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

      <div style={{ width: `${width}px` }}>
        <VegaEmbed
          ref={vegaRef}
          data={state.data}
          spec={initialSpec}
          signalListeners={state.signalListeners}
          dataListeners={state.dataListeners}
          signals={{
            detailHeight:
              configStore.groupKey === GROUP_KEYS.GROUP_SNV ? 0 : 280,
            hoverBar: { group: configStore.hoverGroup },
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
          height={configStore.groupKey === GROUP_KEYS.GROUP_SNV ? 60 : 380}
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

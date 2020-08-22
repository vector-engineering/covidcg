import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import _ from 'underscore';

import VegaEmbed from '../../react_vega/VegaEmbed';
import EmptyPlot from '../Common/EmptyPlot';
import SkeletonElement from '../Common/SkeletonElement';
import LoadingSpinner from '../Common/LoadingSpinner';
import DropdownButton from '../Buttons/DropdownButton';
import { PlotTitle, PlotOptions, OptionSelectContainer } from './Plot.styles';

import initialSpec from '../../vega_specs/cooccurrence_plot.vg.json';
import { ASYNC_STATES } from '../../constants/UI';
import { PLOT_DOWNLOAD_OPTIONS } from '../../constants/download';
import { NORM_MODES } from '../../constants/plotSettings';
import { GROUPS } from '../../constants/groups';

const PlotContainer = styled.div``;

const CooccurrencePlot = observer(({ width }) => {
  const vegaRef = useRef();
  const { configStore, dataStore, plotSettingsStore, UIStore } = useStores();

  const handleDownloadSelect = (option) => {
    // console.log(option);
    // TODO: use the plot options and configStore options to build a more descriptive filename
    //       something like new_lineages_by_day_S_2020-05-03-2020-05-15_NYC.png...
    if (option === PLOT_DOWNLOAD_OPTIONS.DOWNLOAD_DATA) {
      dataStore.downloadSnvCooccurrence();
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

  const onChangeNormMode = (event) =>
    plotSettingsStore.setCooccurrenceNormMode(event.target.value);

  const handleHoverGroup = (...args) => {
    // Don't fire the action if there's no change
    let hoverGroup = args[1] === null ? null : args[1]['group'];
    if (hoverGroup === configStore.hoverGroup) {
      return;
    }

    // Ignore for the 'None' group
    if (hoverGroup === GROUPS.NONE_GROUP) {
      return;
    }

    configStore.updateHoverGroup(hoverGroup);
  };

  const getCooccurrenceData = () => {
    return JSON.parse(JSON.stringify(dataStore.snvCooccurrence));
  };

  const [state, setState] = useState({
    data: {
      //selectedGroups: JSON.parse(JSON.stringify(configStore.selectedGroups)),
      selectedGroups: [],
      cooccurrence_data: getCooccurrenceData(),
    },
    signalListeners: {
      hoverGroup: _.throttle(handleHoverGroup, 100),
    },
    dataListeners: {
      // selectedLocations: handleSelectedLocations,
      // selectedGroups: handleSelectedGroups,
    },
  });

  useEffect(() => {
    // Only update once the SNV data finished processing
    if (UIStore.snvDataState !== ASYNC_STATES.SUCCEEDED) {
      return;
    }

    setState({
      ...state,
      data: {
        //selectedGroups: JSON.parse(JSON.stringify(configStore.selectedGroups)),
        selectedGroups: [],
        cooccurrence_data: getCooccurrenceData(),
      },
    });
  }, [UIStore.snvDataState]);

  if (UIStore.snvDataState === ASYNC_STATES.STARTED) {
    return (
      <div
        style={{
          paddingTop: '12px',
          paddingRight: '24px',
          paddingLeft: '12px',
          paddingBottom: '24px',
        }}
      >
        <SkeletonElement delay={2} height={100}>
          <LoadingSpinner />
        </SkeletonElement>
      </div>
    );
  }

  if (configStore.selectedGroups.length === 0) {
    return (
      <EmptyPlot height={100}>
        <p>
          No SNVs selected. Please select one or more SNVs from the legend,
          frequency plot, or table.
        </p>
      </EmptyPlot>
    );
  }

  // Signals
  let stackOffset, xLabel, xFormat;
  if (plotSettingsStore.cooccurrenceNormMode === NORM_MODES.NORM_COUNTS) {
    stackOffset = 'zero';
    xLabel = 'SNV Frequency';
    xFormat = 's';
  } else if (
    plotSettingsStore.cooccurrenceNormMode === NORM_MODES.NORM_PERCENTAGES
  ) {
    stackOffset = 'normalize';
    xLabel = 'SNV Percentages';
    xFormat = '%';
  }

  // Subtitle text
  const maxShownSnvs = 4;
  let subtitle =
    '(' +
    configStore.selectedGroups
      .slice(0, maxShownSnvs)
      .map((item) => item.group)
      .join(', ');
  if (configStore.selectedGroups.length > maxShownSnvs) {
    subtitle += ', ...)';
  } else {
    subtitle += ')';
  }

  return (
    <PlotContainer>
      <PlotOptions>
        <PlotTitle>
          <span className="title">
            {configStore.getGroupLabel()} Co-occurrence
          </span>
          <span className="subtitle">{subtitle}</span>
        </PlotTitle>
        Show SNVs as{' '}
        <OptionSelectContainer>
          <label>
            <select
              value={plotSettingsStore.cooccurrenceNormMode}
              onChange={onChangeNormMode}
            >
              <option value={NORM_MODES.NORM_COUNTS}>Counts</option>
              <option value={NORM_MODES.NORM_PERCENTAGES}>Percentages</option>
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
      <VegaEmbed
        ref={vegaRef}
        width={width}
        spec={initialSpec}
        data={state.data}
        signalListeners={state.signalListeners}
        dataListeners={state.dataListeners}
        signals={{
          hoverGroup: { group: configStore.hoverGroup },
          stackOffset,
          xLabel,
          xFormat,
        }}
        actions={false}
      />
    </PlotContainer>
  );
});
CooccurrencePlot.propTypes = {
  width: PropTypes.number,
};
CooccurrencePlot.defaultProps = {
  width: 100,
};
export default CooccurrencePlot;

import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { toJS } from 'mobx';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { aggregate } from '../../utils/transform';
import { throttle } from '../../utils/func';

import VegaEmbed from '../../react_vega/VegaEmbed';
import EmptyPlot from '../Common/EmptyPlot';
import SkeletonElement from '../Common/SkeletonElement';
import DropdownButton from '../Buttons/DropdownButton';
import { PlotTitle, PlotOptions, OptionSelectContainer } from './Plot.styles';

import initialSpec from '../../vega_specs/cooccurrence_plot.vg.json';
import {
  ASYNC_STATES,
  PLOT_DOWNLOAD_OPTIONS,
  NORM_MODES,
  GROUPS,
  DNA_OR_AA,
} from '../../constants/defs.json';
import { formatSnv } from '../../utils/snpUtils';

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
    configStore.updateHoverGroup(args[1] === null ? null : args[1]['group']);
  };

  const handleSelectedGroups = (...args) => {
    const curSelectedGroups = toJS(configStore.selectedGroups).map(
      (item) => item.group
    );
    // console.log(curSelectedGroups);

    // Some of the groups might be multiple SNVs, where the
    // group string is the two SNVs separated by ' + '
    // So break those up and then combine into one big list
    // Also make sure we don't process "Reference" or "Other"
    // Array -> Set -> Array to remove duplicates
    const newSelectedGroups = Array.from(
      new Set(
        args[1]
          .map((item) => item.group.split(' + '))
          .reduce((memo, arr) => memo.concat(arr), [])
          .filter((item) => {
            return (
              item !== GROUPS.REFERENCE_GROUP && item !== GROUPS.OTHER_GROUP
            );
          })
      )
    );
    // console.log(newSelectedGroups);

    const groupsInCommon = newSelectedGroups
      .slice()
      .filter((group) => curSelectedGroups.includes(group));

    // The new selection to push to the store is
    // (cur u new) - (cur n new)
    // i.e., the union minus the intersection (XOR)
    // Array -> Set -> Array to remove duplicates
    // Then wrap in an object of { group: group }
    const pushSelectedGroups = Array.from(
      new Set(
        curSelectedGroups
          .concat(newSelectedGroups)
          .filter((group) => !groupsInCommon.includes(group))
      )
    ).map((group) => {
      return { group };
    });
    // console.log(pushSelectedGroups);

    configStore.updateSelectedGroups(pushSelectedGroups);
  };

  const getCooccurrenceData = () => {
    // console.log('GET COOCCURRENCE DATA');
    let newCooccurrenceData = toJS(dataStore.snvCooccurrence);
    newCooccurrenceData = aggregate({
      data: newCooccurrenceData,
      groupby: ['combi', 'snv', 'snvName'],
      fields: ['combiName', 'snvName', 'color', 'count'],
      ops: ['max', 'max', 'max', 'sum'],
      as: ['combiName', 'snvName', 'color', 'count'],
    });

    return newCooccurrenceData;
  };

  const [state, setState] = useState({
    data: {},
    hoverGroup: null,
    signalListeners: {
      hoverGroup: throttle(handleHoverGroup, 50),
    },
    dataListeners: {
      selectedGroups: handleSelectedGroups,
    },
  });

  useEffect(() => {
    setState({
      ...state,
      hoverGroup: { group: configStore.hoverGroup },
    });
  }, [configStore.hoverGroup]);

  const refreshData = () => {
    // Only update once the SNV data finished processing
    if (UIStore.cooccurrenceDataState !== ASYNC_STATES.SUCCEEDED) {
      return;
    }

    setState({
      ...state,
      data: {
        //selectedGroups: toJS(configStore.selectedGroups),
        selectedGroups: [],
        cooccurrence_data: getCooccurrenceData(),
      },
    });
  };

  // Refresh data on mount (i.e., tab change) or when data state changes
  useEffect(refreshData, [UIStore.cooccurrenceDataState]);
  useEffect(refreshData, []);

  if (UIStore.cooccurrenceDataState === ASYNC_STATES.STARTED) {
    return (
      <div
        style={{
          paddingTop: '12px',
          paddingRight: '24px',
          paddingLeft: '12px',
          paddingBottom: '24px',
        }}
      >
        <SkeletonElement delay={2} height={70} />
      </div>
    );
  }

  if (configStore.selectedGroups.length === 0) {
    return (
      <EmptyPlot height={70}>
        <p>
          No SNVs selected. Please select one or more SNVs from the legend,
          frequency plot, or table.
        </p>
      </EmptyPlot>
    );
  } else if (dataStore.snvCooccurrence.length === 0) {
    return (
      <EmptyPlot height={70}>
        <p>
          No SNVs that co-occur with selected SNVs{' '}
          {configStore.selectedGroups
            .map((item) => formatSnv(item.group, configStore.dnaOrAa))
            .join(', ')}
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
    xLabel = 'SNV Frequency (Normalized)';
    xFormat = 's';
  }

  // Subtitle text
  const maxShownSnvs = 4;
  let subtitle = configStore.selectedGroups
    .slice(0, maxShownSnvs)
    .map((item) => formatSnv(item.group, configStore.dnaOrAa))
    .join(', ');
  if (configStore.selectedGroups.length > maxShownSnvs) {
    subtitle += ', ...';
  } else {
    subtitle += '';
  }

  return (
    <PlotContainer>
      <PlotOptions>
        <PlotTitle>
          <span className="title">Co-occurring SNVs of {subtitle}</span>
        </PlotTitle>
        Show SNVs as{' '}
        <OptionSelectContainer>
          <label>
            <select
              value={plotSettingsStore.cooccurrenceNormMode}
              onChange={onChangeNormMode}
            >
              <option value={NORM_MODES.NORM_COUNTS}>Counts</option>
              <option value={NORM_MODES.NORM_PERCENTAGES}>
                Normalized Counts
              </option>
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
          dna: configStore.dnaOrAa === DNA_OR_AA.DNA,
          hoverGroup: state.hoverGroup,
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

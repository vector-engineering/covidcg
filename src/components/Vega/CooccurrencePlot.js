import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { aggregate } from '../../utils/transform';

import _ from 'underscore';

import VegaEmbed from '../../react_vega/VegaEmbed';
import EmptyPlot from '../Common/EmptyPlot';
import SkeletonElement from '../Common/SkeletonElement';
import DropdownButton from '../Buttons/DropdownButton';
import { PlotTitle, PlotOptions, OptionSelectContainer } from './Plot.styles';

import initialSpec from '../../vega_specs/cooccurrence_plot.vg.json';
import { ASYNC_STATES } from '../../constants/UI';
import { PLOT_DOWNLOAD_OPTIONS } from '../../constants/download';
import { NORM_MODES } from '../../constants/plotSettings';
import { GROUPS } from '../../constants/groups';
import { DNA_OR_AA } from '../../constants/config';
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

  const handleHoverGroup = _.debounce((...args) => {
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
  }, 20);

  const handleSelectedGroups = (...args) => {
    const curSelectedGroups = getSelectedGroups().map((item) => item.group);
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
    let newCooccurrenceData = JSON.parse(
      JSON.stringify(dataStore.snvCooccurrence)
    );
    newCooccurrenceData = aggregate({
      data: newCooccurrenceData,
      groupby: ['combi', 'snv', 'snvName'],
      fields: ['combiName', 'snvName', 'color', 'count'],
      ops: ['max', 'max', 'max', 'sum'],
      as: ['combiName', 'snvName', 'color', 'count'],
    });

    return newCooccurrenceData;
  };

  const getSelectedGroups = () => {
    return JSON.parse(JSON.stringify(configStore.selectedGroups));
  };

  const [state, setState] = useState({
    data: {
      //selectedGroups: getSelectedGroups(),
      selectedGroups: [],
      cooccurrence_data: getCooccurrenceData(),
    },
    signalListeners: {
      hoverGroup: _.debounce(handleHoverGroup, 20),
    },
    dataListeners: {
      selectedGroups: handleSelectedGroups,
    },
  });

  useEffect(() => {
    // Only update once the SNV data finished processing
    if (UIStore.cooccurrenceDataState !== ASYNC_STATES.SUCCEEDED) {
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
  }, [UIStore.cooccurrenceDataState]);

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

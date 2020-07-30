import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
// import { toJS } from 'mobx';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import { asyncStates } from '../../stores/uiStore';
import _ from 'underscore';

import EmptyPlot from '../Common/EmptyPlot';
import WarningBox from '../Common/WarningBox';
import DropdownButton from '../Buttons/DropdownButton';
import VegaEmbed from '../../react_vega/VegaEmbed';
import SkeletonElement from '../Common/SkeletonElement';
import LoadingSpinner from '../Common/LoadingSpinner';

// import areaStackSpecInitial from '../vega/area_stack.vl.json';
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

const VegaStackedBars = observer(({ width }) => {
  const vegaRef = useRef();
  const { covidStore, uiStore } = useStores();

  const handleBrush = (...args) => {
    let dateRange = args[1];
    if (dateRange !== null) {
      covidStore.selectDateRange([
        dateRange[0].getTime(),
        dateRange[1].getTime(),
      ]);
    } else {
      // Reset time range
      covidStore.selectDateRange([-1, -1]);
    }
  };

  const handleHoverGroup = (...args) => {
    // Don't fire the action if there's no change
    let hoverGroup = args[1] === null ? null : args[1]['group'];
    if (hoverGroup === covidStore.hoverGroup) {
      return;
    }

    // console.log('Updating store hoverGroup from plot', hoverGroup);
    // console.log('state hovergroup', state.hoverGroup);
    covidStore.updateHoverGroup(hoverGroup);
  };

  const handleSelected = (...args) => {
    // console.log(args);

    // Don't fire if the selection is the same
    if (_.isEqual(args[1], covidStore.selectedGroups)) {
      return;
    }

    covidStore.updateSelectedGroups(args[1]);
  };

  const handleDownloadSelect = (option) => {
    // console.log(option);
    // TODO: use the plot options and covidStore options to build a more descriptive filename
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

  const [state, setState] = useState({
    showWarning: true,
    data: {
      cases_by_date_and_group: JSON.parse(JSON.stringify(covidStore.caseData)),
      selected: JSON.parse(JSON.stringify(covidStore.selectedGroups)),
    },
    signalListeners: {
      detailDomain: _.debounce(handleBrush, 500),
      hoverBar: _.throttle(handleHoverGroup, 100),
    },
    dataListeners: {
      selected: handleSelected,
    },
    areaStackMode: 'counts', // 'percentages' or 'counts'
    countMode: 'new', // 'new' or 'cumulative'
    dateBin: 'day', // 'day', 'week', 'month'
  });

  const onDismissWarning = () => {
    setState({
      ...state,
      showWarning: false,
    });
  };

  const onChangeAreaStackMode = (event) =>
    setState({ ...state, areaStackMode: event.target.value });
  const onChangeCountMode = (event) =>
    setState({ ...state, countMode: event.target.value });
  const onChangeDateBin = (event) =>
    setState({ ...state, dateBin: event.target.value });

  // Update internal caseData copy
  useEffect(() => {
    // console.log('Update caseData');

    // let newCaseData;

    // if (covidStore.groupsToKeep) {
    //   newCaseData = [];
    //   const addedOthersByDate = {};

    //   covidStore.caseData.forEach((row, key) => {
    //     if (!covidStore.groupsToKeep[row.group]) {
    //       row.group = 'other';
    //       row.color = '#cccccc';
    //       if (addedOthersByDate[row.date]) {
    //         covidStore.caseData[addedOthersByDate[row.date]].cases_sum +=
    //           row.cases_sum;
    //       } else {
    //         addedOthersByDate[row.date] = key;
    //         newCaseData.push(row);
    //       }
    //     } else {
    //       newCaseData.push(row);
    //     }
    //   });
    // }

    setState({
      ...state,
      data: {
        ...state.data,
        cases_by_date_and_group: JSON.parse(
          JSON.stringify(covidStore.caseData)
        ),
      },
    });
  }, [covidStore.caseData]);

  // Update internal selected groups copy
  useEffect(() => {
    // console.log('Update selectedGroups');
    setState({
      ...state,
      data: {
        ...state.data,
        selected: JSON.parse(JSON.stringify(covidStore.selectedGroups)),
      },
    });
  }, [covidStore.selectedGroups]);

  // For development in Vega Editor
  // console.log(JSON.stringify(caseData));

  let areaStackTitle = '';
  if (state.countMode === 'cumulative') {
    areaStackTitle += 'Cumulative ';
  } else if (state.countMode === 'new') {
    areaStackTitle += 'New ';
  }

  if (covidStore.groupKey === 'lineage') {
    areaStackTitle += 'Lineage ';
  } else if (covidStore.groupKey === 'clade') {
    areaStackTitle += 'Clade ';
  } else if (covidStore.groupKey === 'snp') {
    if (covidStore.dnaOrAa === 'dna') {
      areaStackTitle += 'NT';
    } else {
      areaStackTitle += 'AA';
    }
    areaStackTitle += ' SNV ';
  }
  areaStackTitle +=
    state.areaStackMode === 'percentages' ? 'Percentages' : 'Counts';

  if (state.dateBin === 'day') {
    areaStackTitle += ' by Day';
  } else if (state.dateBin === 'week') {
    areaStackTitle += ' by Week';
  } else if (state.dateBin === 'month') {
    areaStackTitle += ' by Month';
  }

  // Set the stack offset mode
  const stackOffset =
    state.areaStackMode === 'percentages' ? 'normalize' : 'zero';
  // Set the date bin
  let dateBin;
  if (state.dateBin === 'day') {
    dateBin = 1000 * 60 * 60 * 24;
  } else if (state.dateBin === 'week') {
    dateBin = 1000 * 60 * 60 * 24 * 7;
  } else if (state.dateBin === 'month') {
    dateBin = 1000 * 60 * 60 * 24 * 30;
  }
  // If running in cumulative mode, add the vega transformation
  // By default the cumulative transformation is dumped into a column
  // "cases_sum_cumulative", so if active, just overwrite the "cases_sum"
  // column with this cumulative count
  const cumulativeWindow =
    state.countMode === 'cumulative' ? [null, 0] : [0, 0];

  // Adapt labels to groupings
  let detailYLabel = '';
  if (state.countMode === 'cumulative') {
    detailYLabel += 'Cumulative ';
  }
  if (state.areaStackMode === 'percentages') {
    detailYLabel += '% ';
  }
  if (covidStore.groupKey === 'lineage') {
    detailYLabel += 'Sequences by Lineage';
  } else if (covidStore.groupKey === 'clade') {
    detailYLabel += 'Sequences by Clade';
  } else if (covidStore.groupKey === 'snp') {
    detailYLabel +=
      'Sequences by ' + (covidStore.dnaOrAa === 'dna' ? 'NT' : 'AA') + ' SNV';
  }

  if (uiStore.caseDataState === asyncStates.STARTED) {
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

  if (covidStore.selectedLocationIds.length === 0) {
    return (
      <EmptyPlot height={250}>
        <p>
          No locations selected. Please select one or more locations from the
          sidebar, under &quot;Selected Locations&quot;, to compare counts of{' '}
          <b>{covidStore.getGroupLabel()}</b> between them.
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
            <select value={state.countMode} onChange={onChangeCountMode}>
              <option value="new">New</option>
              <option value="cumulative">Cumulative</option>
            </select>
          </label>
        </SelectContainer>
        sequences, shown as{' '}
        <SelectContainer>
          <label>
            <select
              value={state.areaStackMode}
              onChange={onChangeAreaStackMode}
            >
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
            hoverBar: { group: covidStore.hoverGroup },
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

VegaStackedBars.propTypes = {
  width: PropTypes.number.isRequired,
};

export default VegaStackedBars;

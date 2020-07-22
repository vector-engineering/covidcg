import React, { useState } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
// import { toJS } from 'mobx';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import _ from 'underscore';

import VegaEmbed from '../../react_vega/VegaEmbed';

// import areaStackSpecInitial from '../vega/area_stack.vl.json';
import barStackSpecInitial from '../../vega/bar_stack_v1.vg.json';

const PlotOptions = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-start;
  margin-bottom: 10px;

  .area-stack-title {
    font-size: 1.25em;
    margin-right: 10px;
    padding-right: 10px;
    padding-left: 18px;

    border-right: 1px solid #ccc;
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
  const { covidStore } = useStores();

  const [state, setState] = useState({
    areaStackMode: 'counts', // 'percentages' or 'counts'
    countMode: 'new', // 'new' or 'cumulative'
    dateBin: 'day', // 'day', 'week', 'month'
  });

  const onChangeAreaStackMode = (event) =>
    setState({ ...state, areaStackMode: event.target.value });
  const onChangeCountMode = (event) =>
    setState({ ...state, countMode: event.target.value });
  const onChangeDateBin = (event) =>
    setState({ ...state, dateBin: event.target.value });

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
    covidStore.updateHoverGroup(hoverGroup);
  };

  const handleSelected = (...args) => {
    // console.log(args);
    covidStore.updateSelectedGroups(args[1]);
  };

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

  // Make a deep copy of the Vega spec so we can edit it
  const barStackSpec = JSON.parse(JSON.stringify(barStackSpecInitial));

  // Set the width manually
  barStackSpec['width'] = width;
  barStackSpec['marks'][0]['encode']['enter']['width']['value'] = width;
  barStackSpec['marks'][1]['encode']['enter']['width']['value'] = width;

  // TODO: these are signals and should be able to be set when passed
  //       through the signal prop object. but for some reason it doesn't
  //       trigger the proper re-render, probably because I have no idea
  //       how the Vega View API actually works. For now manually modifying
  //       the spec and triggering a hard re-render will work...
  // Set the stack offset mode
  let stackOffsetSignal = _.findWhere(barStackSpec['signals'], {
    name: 'stackOffset',
  });
  if (state.areaStackMode === 'percentages') {
    stackOffsetSignal['value'] = 'normalize';
  } else {
    stackOffsetSignal['value'] = 'zero';
  }

  // Set the date bin
  let dateBinSignal = _.findWhere(barStackSpec['signals'], {
    name: 'dateBin',
  });
  if (state.dateBin === 'day') {
    dateBinSignal['value'] = 1000 * 60 * 60 * 24;
  } else if (state.dateBin === 'week') {
    dateBinSignal['value'] = 1000 * 60 * 60 * 24 * 7;
  } else if (state.dateBin === 'month') {
    dateBinSignal['value'] = 1000 * 60 * 60 * 24 * 30;
  }

  // If running in cumulative mode, add the vega transformation
  // By default the cumulative transformation is dumped into a column
  // "cases_sum_cumulative", so if active, just overwrite the "cases_sum"
  // column with this cumulative count
  let windowFieldSignal = _.findWhere(barStackSpec['signals'], {
    name: 'windowField',
  });
  if (state.countMode === 'cumulative') {
    windowFieldSignal['value'] = 'cases_sum';
  } else {
    windowFieldSignal['value'] = 'cases_sum_cumulative';
  }

  // Adapt labels to groupings
  let detailYLabelSignal = _.findWhere(barStackSpec['signals'], {
    name: 'detailYLabel',
  });
  // let detailBars = detailGroup['marks'][0]['marks'][0];
  let detailYLabel = '';
  if (state.countMode === 'cumulative') {
    detailYLabel += 'Cumulative ';
  }
  if (state.areaStackMode === 'percentages') {
    detailYLabel += '% ';
  }
  if (covidStore.groupKey === 'lineage') {
    detailYLabel += 'Sequences by Lineage';
  } else if (covidStore.groupKey === 'snp') {
    detailYLabel +=
      'Sequences by ' + (covidStore.dnaOrAa === 'dna' ? 'NT' : 'AA') + ' SNP';
  }
  detailYLabelSignal['value'] = detailYLabel;

  let caseData = JSON.parse(JSON.stringify(covidStore.caseData));
  let selectedGroups = JSON.parse(JSON.stringify(covidStore.selectedGroups));

  // Compute counts for each group
  let countsPerGroup = {};
  caseData.forEach((row) => {
    if (!Object.prototype.hasOwnProperty.call(countsPerGroup, row.group)) {
      countsPerGroup[row.group] = 0;
    }
    countsPerGroup[row.group] += row.cases_sum;
  });

  // Map counts per group back onto main dataset
  caseData.forEach((row) => {
    row['group_counts'] = countsPerGroup[row.group];
  });

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
  } else if (covidStore.groupKey === 'snp') {
    if (covidStore.dnaOrAa === 'dna') {
      areaStackTitle += 'NT';
    } else {
      areaStackTitle += 'AA';
    }
    areaStackTitle += ' SNP ';
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

  return (
    <div>
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
      </PlotOptions>

      <div style={{ width: `${width - 150}px` }}>
        <VegaEmbed
          data={{
            cases_by_date_and_group: caseData,
            selected: selectedGroups,
          }}
          spec={barStackSpec}
          signalListeners={{
            detailDomain: _.debounce(handleBrush, 500),
            hoverBar: _.throttle(handleHoverGroup, 100),
          }}
          dataListeners={{
            selected: handleSelected,
          }}
          signals={{
            hoverBar: covidStore.hoverGroup,
          }}
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

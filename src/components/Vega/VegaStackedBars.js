import React, { useState, useEffect, useRef } from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';
// import { toJS } from 'mobx';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import _ from 'underscore';

import DropdownButton from '../Buttons/DropdownButton';
import VegaEmbed from '../../react_vega/VegaEmbed';

// import areaStackSpecInitial from '../vega/area_stack.vl.json';
import barStackSpecInitial from '../../vega/bar_stack_v1.vg.json';

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
  const { covidStore } = useStores();

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
    data: {
      cases_by_date_and_group: JSON.parse(JSON.stringify(covidStore.caseData)),
      selected: JSON.parse(JSON.stringify(covidStore.selectedGroups)),
    },
    spec: JSON.parse(JSON.stringify(barStackSpecInitial)),
    signalListeners: {
      detailDomain: _.debounce(handleBrush, 500),
      hoverBar: _.throttle(handleHoverGroup, 100),
    },
    dataListeners: {
      selected: handleSelected,
    },
    hoverGroup: {},
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

  useEffect(() => {
    let spec = JSON.parse(JSON.stringify(state.spec));

    // TODO: these are signals and should be able to be set when passed
    //       through the signal prop object. but for some reason it doesn't
    //       trigger the proper re-render, probably because I have no idea
    //       how the Vega View API actually works. For now manually modifying
    //       the spec and triggering a hard re-render will work...
    // Set the stack offset mode
    let stackOffsetSignal = _.findWhere(spec['signals'], {
      name: 'stackOffset',
    });
    if (state.areaStackMode === 'percentages') {
      stackOffsetSignal['value'] = 'normalize';
    } else {
      stackOffsetSignal['value'] = 'zero';
    }

    // Set the date bin
    let dateBinSignal = _.findWhere(spec['signals'], {
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
    let windowFieldSignal = _.findWhere(spec['signals'], {
      name: 'windowField',
    });
    if (state.countMode === 'cumulative') {
      windowFieldSignal['value'] = 'cases_sum';
    } else {
      windowFieldSignal['value'] = 'cases_sum_cumulative';
    }

    // Adapt labels to groupings
    let detailYLabelSignal = _.findWhere(spec['signals'], {
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
    } else if (covidStore.groupKey === 'clade') {
      detailYLabel += 'Sequences by Clade';
    } else if (covidStore.groupKey === 'snp') {
      detailYLabel +=
        'Sequences by ' + (covidStore.dnaOrAa === 'dna' ? 'NT' : 'AA') + ' SNP';
    }
    detailYLabelSignal['value'] = detailYLabel;

    setState({ ...state, spec });
  }, [
    state.areaStackMode,
    state.dateBin,
    state.countMode,
    covidStore.groupKey,
    covidStore.dnaOrAa,
  ]);

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
          spec={state.spec}
          signalListeners={state.signalListeners}
          dataListeners={state.dataListeners}
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

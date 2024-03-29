import React, { useState, useEffect } from 'react';
import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';

import ReactTooltip from 'react-tooltip';
import QuestionButton from '../Buttons/QuestionButton';

import { LOW_FREQ_FILTER_TYPES } from '../../constants/defs.json';

import {
  LowFreqFilterContainer,
  LowFreqFilterSelectContainer,
  LowFreqFilterInputContainer,
  ApplyButton,
} from './LowFreqFilter.styles';

const LowFreqFilter = observer(
  ({
    lowFreqFilterType,
    lowFreqFilterValue,
    updateLowFreqFilterType,
    updateLowFreqFilterValue,
    ...other
  }) => {
    const { configStore } = useStores();
    const [state, setState] = useState({
      filterType: lowFreqFilterType,
      filterValue: lowFreqFilterValue,
      hasChanged: false,
    });

    const prefix = {
      [LOW_FREQ_FILTER_TYPES.GROUP_COUNTS]: 'Top',
      [LOW_FREQ_FILTER_TYPES.LOCAL_COUNTS]: 'Minimum',
    };
    const suffix = {
      [LOW_FREQ_FILTER_TYPES.GROUP_COUNTS]: configStore.getGroupLabel() + 's',
      [LOW_FREQ_FILTER_TYPES.LOCAL_COUNTS]: 'counts',
    };
    const defaults = {
      [LOW_FREQ_FILTER_TYPES.GROUP_COUNTS]: 20,
      [LOW_FREQ_FILTER_TYPES.LOCAL_COUNTS]: 50,
    };

    const onChangeLowFreqFilter = (e) => {
      const newFilterType = e.target.value;
      if (newFilterType === state.filterType) {
        return;
      }

      setState({
        filterType: newFilterType,
        filterValue: defaults[newFilterType],
        hasChanged: true,
      });
    };
    const onChangeLowFreqValue = (e) => {
      const newFilterValue = e.target.value;
      setState({
        ...state,
        filterValue: newFilterValue,
        hasChanged: !(newFilterValue === state.filterValue),
      });
    };

    const flushChanges = () => {
      updateLowFreqFilterType(state.filterType);
      updateLowFreqFilterValue(state.filterValue);
      setState({
        ...state,
        hasChanged: false,
      });
    };

    useEffect(() => {
      ReactTooltip.rebuild();
    }, []);

    return (
      <LowFreqFilterContainer>
        <LowFreqFilterSelectContainer>
          Collapse {configStore.getGroupLabel()}s by:
          <select value={state.filterType} onChange={onChangeLowFreqFilter}>
            <option value={LOW_FREQ_FILTER_TYPES.GROUP_COUNTS}>Top N</option>
            <option value={LOW_FREQ_FILTER_TYPES.LOCAL_COUNTS}>
              Minimum Counts
            </option>
          </select>
        </LowFreqFilterSelectContainer>
        <LowFreqFilterInputContainer>
          {prefix[state.filterType] + ' '}
          <input
            value={state.filterValue}
            onChange={onChangeLowFreqValue}
            {...other}
          ></input>
          {' ' + suffix[state.filterType]}
        </LowFreqFilterInputContainer>
        <ApplyButton
          disabled={!state.hasChanged}
          invalid={!state.hasChanged}
          onClick={flushChanges}
        >
          Apply
        </ApplyButton>
        <QuestionButton
          data-tip={`
            <p>${configStore.getGroupLabel(
              configStore.groupKey,
              configStore.dnaOrAa
            )}s that do not meet the following criteria will be grouped into the "Other" group. This is done to increase performance in the app</p><p>Including more groups gives more detail into the data, but may come at the cost of app performance.</p>
            <p>For the "Minimum Counts" option, <b>${configStore.getGroupLabel(
              configStore.groupKey,
              configStore.dnaOrAa
            )}s</b> with less than N counts in the selected locations will be grouped into "Other"</p>
            <p>For the "Top N" option, only show the top N <b>${configStore.getGroupLabel(
              configStore.groupKey,
              configStore.dnaOrAa
            )}s</b> by counts in the selected locations</p>
          `}
          data-html="true"
          data-for="main-tooltip"
        />
      </LowFreqFilterContainer>
    );
  }
);

LowFreqFilter.propTypes = {
  lowFreqFilterType: PropTypes.oneOf([
    LOW_FREQ_FILTER_TYPES.GROUP_COUNTS,
    LOW_FREQ_FILTER_TYPES.LOCAL_COUNTS,
  ]).isRequired,
  lowFreqFilterValue: PropTypes.oneOfType([PropTypes.number, PropTypes.string])
    .isRequired,
  updateLowFreqFilterType: PropTypes.func.isRequired,
  updateLowFreqFilterValue: PropTypes.func.isRequired,
};

export default LowFreqFilter;

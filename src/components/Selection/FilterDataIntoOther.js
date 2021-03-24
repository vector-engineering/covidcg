import React from 'react';

import QuestionButton from '../Buttons/QuestionButton';

import {
  SelectContainer,
  SelectItem,
  Radio,
  RadioLabel,
  NumberInput,
} from './FilterDataIntoOther.styles';

import { useStores } from '../../stores/connect';
import { observer } from 'mobx-react';
import { LOW_FREQ_FILTER_TYPES } from '../../constants/defs.json';

const FilterDataIntoOther = observer(
  ({
    groupKey,
    dnaOrAa,

    lowFreqFilterType,
    minLocalCounts,
    minGlobalCounts,
    maxGroupCounts,

    setLowFreqFilterType,
    setMinLocalCounts,
    setMinGlobalCounts,
    setMaxGroupCounts,
  }) => {
    const { configStore } = useStores();

    const _setLowFreqFilterType = (e) => {
      setLowFreqFilterType(e.target.value);
    };

    const _setMinLocalCounts = (e) => {
      setMinLocalCounts(e.target.value);
    };

    const _setMinGlobalCounts = (e) => {
      setMinGlobalCounts(e.target.value);
    };

    const _setMaxGroupCounts = (e) => {
      setMaxGroupCounts(e.target.value);
    };

    return (
      <SelectContainer>
        <span className="title">
          Filter Low Frequency {configStore.getGroupLabel(groupKey, dnaOrAa)}s
          <QuestionButton
            data-tip={`<p>${configStore.getGroupLabel(
              groupKey,
              dnaOrAa
            )}s that do not meet the following criteria will be grouped into the "Other" group. This is done to increase performance in the app</p><p>Including more groups gives more detail into the data, but may come at the cost of app performance.</p>`}
            data-html={true}
            data-for="main-tooltip"
          />
        </span>
        <SelectItem>
          <Radio>
            <input
              type="radio"
              value={LOW_FREQ_FILTER_TYPES.LOCAL_COUNTS}
              onChange={_setLowFreqFilterType}
              checked={lowFreqFilterType === LOW_FREQ_FILTER_TYPES.LOCAL_COUNTS}
            />
            <RadioLabel
              itemSelected={
                lowFreqFilterType === LOW_FREQ_FILTER_TYPES.LOCAL_COUNTS
              }
            >
              Minimum local counts
            </RadioLabel>
          </Radio>
          <NumberInput>
            <input
              type="number"
              value={minLocalCounts}
              onChange={_setMinLocalCounts}
              disabled={
                lowFreqFilterType !== LOW_FREQ_FILTER_TYPES.LOCAL_COUNTS
              }
              min={1}
              step={1}
            />
          </NumberInput>
          <QuestionButton
            data-tip={`<p><b>${configStore.getGroupLabel(
              groupKey,
              dnaOrAa
            )}s</b> with less than <b>${minLocalCounts}</b> counts in the selected locations will be grouped into "Other"</p>`}
            data-html="true"
            data-for="main-tooltip"
          />
        </SelectItem>
        <SelectItem>
          <Radio>
            <input
              type="radio"
              value={LOW_FREQ_FILTER_TYPES.GLOBAL_COUNTS}
              onChange={_setLowFreqFilterType}
              checked={
                lowFreqFilterType === LOW_FREQ_FILTER_TYPES.GLOBAL_COUNTS
              }
            />
            <RadioLabel
              itemSelected={
                lowFreqFilterType === LOW_FREQ_FILTER_TYPES.GLOBAL_COUNTS
              }
            >
              Minimum global counts
            </RadioLabel>
          </Radio>
          <NumberInput>
            <input
              type="number"
              value={minGlobalCounts}
              onChange={_setMinGlobalCounts}
              disabled={
                lowFreqFilterType !== LOW_FREQ_FILTER_TYPES.GLOBAL_COUNTS
              }
              min={1}
              step={1}
            />
          </NumberInput>
          <QuestionButton
            data-tip={`<p><b>${configStore.getGroupLabel(
              groupKey,
              dnaOrAa
            )}s</b> with less than <b>${minGlobalCounts}</b> counts globally will be grouped into "Other"</p>`}
            data-html="true"
            data-for="main-tooltip"
          />
        </SelectItem>
        <SelectItem>
          <Radio>
            <input
              type="radio"
              value={LOW_FREQ_FILTER_TYPES.GROUP_COUNTS}
              onChange={_setLowFreqFilterType}
              checked={lowFreqFilterType === LOW_FREQ_FILTER_TYPES.GROUP_COUNTS}
            />
            <RadioLabel
              itemSelected={
                lowFreqFilterType === LOW_FREQ_FILTER_TYPES.GROUP_COUNTS
              }
            >
              Show Top N {configStore.getGroupLabel(groupKey, dnaOrAa)}s
            </RadioLabel>
          </Radio>
          <NumberInput>
            <input
              type="number"
              value={maxGroupCounts}
              onChange={_setMaxGroupCounts}
              disabled={
                lowFreqFilterType !== LOW_FREQ_FILTER_TYPES.GROUP_COUNTS
              }
              min={1}
              step={1}
            />
          </NumberInput>
          <QuestionButton
            data-tip={`<p>Only show the top <b>${maxGroupCounts}</b> <b>${configStore.getGroupLabel(
              groupKey,
              dnaOrAa
            )}s</b> by counts in the selected locations</p>`}
            data-html="true"
            data-for="main-tooltip"
          />
        </SelectItem>
      </SelectContainer>
    );
  }
);

export default FilterDataIntoOther;

import React from 'react';

import ReactTooltip from 'react-tooltip';
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

const FilterDataIntoOther = observer(() => {
  const { configStore } = useStores();

  const setFilterType = (e) => {
    configStore.setLowFreqFilterType(e.target.value);
  };

  const setMinLocalCountsToShow = (e) => {
    configStore.setMinLocalCountsToShow(e.target.value);
  };

  const setMinGlobalCountsToShow = (e) => {
    configStore.setMinGlobalCountsToShow(e.target.value);
  };

  const setMaxGroupCounts = (e) => {
    configStore.setMaxGroupCounts(e.target.value);
  };

  return (
    <SelectContainer>
      <ReactTooltip
        id="low-count-filter-tooltip"
        type="light"
        effect="solid"
        border={true}
        borderColor="#888"
      />
      <span className="title">
        Filter Low Frequency {configStore.getGroupLabel()}s
        <QuestionButton
          data-tip={`<p>${configStore.getGroupLabel()}s that do not meet the following criteria will be grouped into the "Other" group. This is done to increase performance in the app</p><p>Including more groups gives more detail into the data, but may come at the cost of app performance.</p>`}
          data-html={true}
          data-for="low-count-filter-tooltip"
        />
      </span>
      <SelectItem>
        <Radio>
          <input
            type="radio"
            value={LOW_FREQ_FILTER_TYPES.LOCAL_COUNTS}
            onChange={setFilterType}
            checked={
              configStore.lowFreqFilterType ===
              LOW_FREQ_FILTER_TYPES.LOCAL_COUNTS
            }
          />
          <RadioLabel
            itemSelected={
              configStore.lowFreqFilterType ===
              LOW_FREQ_FILTER_TYPES.LOCAL_COUNTS
            }
          >
            Minimum local counts
          </RadioLabel>
        </Radio>
        <NumberInput>
          <input
            type="number"
            value={configStore.minLocalCountsToShow}
            onChange={setMinLocalCountsToShow}
            disabled={
              configStore.lowFreqFilterType !==
              LOW_FREQ_FILTER_TYPES.LOCAL_COUNTS
            }
            min={1}
            step={1}
          />
        </NumberInput>
        <QuestionButton
          data-tip={`<p><b>${configStore.getGroupLabel()}s</b> with less than <b>${
            configStore.minLocalCountsToShow
          }</b> counts in the selected locations will be grouped into "Other"</p>`}
          data-html="true"
          data-for="low-count-filter-tooltip"
        />
      </SelectItem>
      <SelectItem>
        <Radio>
          <input
            type="radio"
            value={LOW_FREQ_FILTER_TYPES.GLOBAL_COUNTS}
            onChange={setFilterType}
            checked={
              configStore.lowFreqFilterType ===
              LOW_FREQ_FILTER_TYPES.GLOBAL_COUNTS
            }
          />
          <RadioLabel
            itemSelected={
              configStore.lowFreqFilterType ===
              LOW_FREQ_FILTER_TYPES.GLOBAL_COUNTS
            }
          >
            Minimum global counts
          </RadioLabel>
        </Radio>
        <NumberInput>
          <input
            type="number"
            value={configStore.minGlobalCountsToShow}
            onChange={setMinGlobalCountsToShow}
            disabled={
              configStore.lowFreqFilterType !==
              LOW_FREQ_FILTER_TYPES.GLOBAL_COUNTS
            }
            min={1}
            step={1}
          />
        </NumberInput>
        <QuestionButton
          data-tip={`<p><b>${configStore.getGroupLabel()}s</b> with less than <b>${
            configStore.minGlobalCountsToShow
          }</b> counts globally will be grouped into "Other"</p>`}
          data-html="true"
          data-for="low-count-filter-tooltip"
        />
      </SelectItem>
      <SelectItem>
        <Radio>
          <input
            type="radio"
            value={LOW_FREQ_FILTER_TYPES.GROUP_COUNTS}
            onChange={setFilterType}
            checked={
              configStore.lowFreqFilterType ===
              LOW_FREQ_FILTER_TYPES.GROUP_COUNTS
            }
          />
          <RadioLabel
            itemSelected={
              configStore.lowFreqFilterType ===
              LOW_FREQ_FILTER_TYPES.GROUP_COUNTS
            }
          >
            Show Top N {configStore.getGroupLabel()}s
          </RadioLabel>
        </Radio>
        <NumberInput>
          <input
            type="number"
            value={configStore.maxGroupCounts}
            onChange={setMaxGroupCounts}
            disabled={
              configStore.lowFreqFilterType !==
              LOW_FREQ_FILTER_TYPES.GROUP_COUNTS
            }
            min={1}
            step={1}
          />
        </NumberInput>
        <QuestionButton
          data-tip={`<p>Only show the top <b>${
            configStore.maxGroupCounts
          }</b> <b>${configStore.getGroupLabel()}s</b> by counts in the selected locations</p>`}
          data-html="true"
          data-for="low-count-filter-tooltip"
        />
      </SelectItem>
    </SelectContainer>
  );
});

export default FilterDataIntoOther;

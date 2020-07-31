import React from 'react';
import styled from 'styled-components';

import QuestionButton from '../Buttons/QuestionButton';

import { useStores } from '../../stores/connect';
import { observer } from 'mobx-react';
import { LOW_FREQ_FILTER_TYPES } from '../../stores/configStore';

export const SelectContainer = styled.div`
  padding: 3px 13px 5px 13px;

  display: flex;
  flex-direction: column;
`;

const SelectItem = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;

  margin-bottom: 3px;
`;

const Radio = styled.label`
  flex-shrink: 0;
  margin-right: 5px;
`;

const Spacer = styled.div`
  flex-grow: 1;
`;

const NumberInput = styled.label`
  flex-shrink: 0;
  width: 50px;
  margin-right: 5px;

  input {
    width: 100%;
  }
`;

const FilterDataIntoOther = observer(() => {
  const { configStore } = useStores();

  const setFilterType = (e) => {
    configStore.setLowFreqFilterType(e.target.value);
  };

  return (
    <SelectContainer>
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
          Minimum local counts
        </Radio>
        <Spacer />
        <NumberInput>
          <input
            type="number"
            value={configStore.minLocalCountsToShow}
            onChange={(e) => configStore.setMinLocalCounts(e.target.value)}
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
          data-for="tooltip-filter-sidebar"
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
          Minimum global counts
        </Radio>
        <Spacer />
        <NumberInput>
          <input
            type="number"
            value={configStore.minGlobalCountsToShow}
            onChange={(e) => configStore.setMinGlobalCounts(e.target.value)}
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
          data-for="tooltip-filter-sidebar"
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
          Show Top N {configStore.getGroupLabel().toLowerCase()}s
        </Radio>
        <Spacer />
        <NumberInput>
          <input
            type="number"
            value={configStore.maxLineagesToShow}
            onChange={(e) => configStore.setMaxLineages(e.target.value)}
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
            configStore.maxLineagesToShow
          }</b> <b>${configStore.getGroupLabel()}s</b> by counts in the selected locations</p>`}
          data-html="true"
          data-for="tooltip-filter-sidebar"
        />
      </SelectItem>
    </SelectContainer>
  );
});

export default FilterDataIntoOther;

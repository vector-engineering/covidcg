import React from 'react';
import styled from 'styled-components';

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
  width: 80px;

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
          Max shown {configStore.getGroupLabel().toLowerCase()}s
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
      </SelectItem>
    </SelectContainer>
  );
});

export default FilterDataIntoOther;

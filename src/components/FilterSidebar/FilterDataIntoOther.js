import React, { useState } from 'react';
import styled from 'styled-components';

import { useStores } from '../../stores/connect';
import { observer } from 'mobx-react';

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
  const [filterType, setFilterType] = useState('LOCAL_COUNTS');

  return (
    <SelectContainer>
      <SelectItem>
        <Radio>
          <input
            type="radio"
            value="LOCAL_COUNTS"
            onChange={(e) => setFilterType(e.target.value)}
            checked={filterType === 'LOCAL_COUNTS'}
          />
          Minimum local counts
        </Radio>
        <Spacer />
        <NumberInput>
          <input
            type="number"
            value={configStore.minLocalCountsToShow}
            onChange={(e) => configStore.setMinLocalCounts(e.target.value)}
            disabled={filterType !== 'LOCAL_COUNTS'}
            min={1}
            step={1}
          />
        </NumberInput>
      </SelectItem>
      <SelectItem>
        <Radio>
          <input
            type="radio"
            value="GROUP_COUNTS"
            onChange={(e) => setFilterType(e.target.value)}
            checked={filterType === 'GROUP_COUNTS'}
          />
          Max shown {configStore.getGroupLabel().toLowerCase()}s
        </Radio>
        <Spacer />
        <NumberInput>
          <input
            type="number"
            value={configStore.maxLineagesToShow}
            onChange={(e) => configStore.setMaxLineages(e.target.value)}
            disabled={filterType !== 'GROUP_COUNTS'}
            min={1}
            step={1}
          />
        </NumberInput>
      </SelectItem>
    </SelectContainer>
  );
});

export default FilterDataIntoOther;

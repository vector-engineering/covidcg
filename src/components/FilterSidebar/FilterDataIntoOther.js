import React, { useState } from 'react';
import styled from 'styled-components';
import Slider from 'react-input-slider';

import { Title } from './SidebarAccordionWrapper';
import { useStores } from '../../stores/connect';
import { observer } from 'mobx-react';

export const SelectContainer = styled.div`
  padding: 7px 13px 5px 13px;
  border-top: 1px solid #aaa;

  display: flex;
  flex-direction: column;
`;

const SliderContainer = styled.div`
  display: flex;
  flex-direction: column;
  padding-left: 24px;
`;

const FilterDataIntoOther = observer(() => {
  const { covidStore } = useStores();
  const [filterType, setFilterType] = useState('LOCAL_COUNTS');

  return (
    <SelectContainer>
      <Title>Collapse Low Freq. Data</Title>
      <label>
        <input
          type="radio"
          value="LOCAL_COUNTS"
          onChange={(e) => setFilterType(e.target.value)}
          checked={filterType === 'LOCAL_COUNTS'}
        />
        Limit lineages by total sequences
      </label>
      {filterType === 'LOCAL_COUNTS' && (
        <SliderContainer>
          {covidStore.minLocalCountsToShow}
          <Slider
            axis="x"
            x={covidStore.minLocalCountsToShow}
            onChange={({ x }) => {
              covidStore.setMinLocalCounts(x);
            }}
            xmin={1}
            xmax={500}
          />
        </SliderContainer>
      )}
      <label>
        <input
          type="radio"
          value="GROUP_COUNTS"
          onChange={(e) => setFilterType(e.target.value)}
          checked={filterType === 'GROUP_COUNTS'}
        />
        Limit number of shown lineages
      </label>
      {filterType === 'GROUP_COUNTS' && (
        <SliderContainer>
          {covidStore.maxLineagesToShow}
          <Slider
            axis="x"
            x={covidStore.maxLineagesToShow}
            onChange={({ x }) => {
              covidStore.setMaxLineages(x);
            }}
            xmin={2}
            xmax={50}
          />
        </SliderContainer>
      )}
    </SelectContainer>
  );
});

export default FilterDataIntoOther;

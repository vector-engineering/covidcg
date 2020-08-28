import React, { useState, useEffect } from 'react';
import styled from 'styled-components';

import Button from '../Buttons/Button';
import QuestionButton from '../Buttons/QuestionButton';

import { useStores } from '../../stores/connect';
import { observer } from 'mobx-react';
import { LOW_FREQ_FILTER_TYPES } from '../../constants/config';

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
  display: flex;
  flex-direction: row;
  align-items: center;

  flex-shrink: 0;
  margin-right: 5px;
`;

const RadioLabel = styled.span`
  color: ${({ itemSelected }) => (itemSelected ? 'inherit' : '#AAA')};
  margin-left: 3px;
`;
RadioLabel.defaultProps = {
  itemSelected: false,
};

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

const ConfirmButton = styled(Button)`
  width: 120px;
  font-size: 14px;
  display: ${({ show }) => (show ? 'block' : 'none')};
`;
ConfirmButton.defaultProps = {
  show: true,
};

const FilterDataIntoOther = observer(() => {
  const { configStore } = useStores();
  const [state, setState] = useState({
    lowFreqFilterType: configStore.lowFreqFilterType,
    minLocalCountsToShow: configStore.minLocalCountsToShow,
    minGlobalCountsToShow: configStore.minGlobalCountsToShow,
    maxGroupCounts: configStore.maxGroupCounts,
    hasChanged: false,
  });

  const setFilterType = (e) => {
    // configStore.setLowFreqFilterType(e.target.value);
    setState({
      ...state,
      localFreqFilterType: e.target.value,
      hasChanged: e.target.value !== configStore.lowFreqFilterType,
    });
  };

  const setMinLocalCountsToShow = (e) => {
    setState({
      ...state,
      minLocalCountsToShow: e.target.value,
      hasChanged: e.target.value !== configStore.minLocalCountsToShow,
    });
  };

  const setMinGlobalCountsToShow = (e) => {
    setState({
      ...state,
      minGlobalCountsToShow: e.target.value,
      hasChanged: e.target.value !== configStore.minGlobalCountsToShow,
    });
  };

  const setMaxGroupCounts = (e) => {
    setState({
      ...state,
      maxGroupCounts: e.target.value,
      hasChanged: e.target.value !== configStore.maxGroupCounts,
    });
  };

  const updateStore = () => {
    configStore.setLowFreqFilters({
      lowFreqFilterType: state.lowFreqFilterType,
      minLocalCountsToShow: state.minLocalCountsToShow,
      minGlobalCountsToShow: state.minGlobalCountsToShow,
      maxGroupCounts: state.maxGroupCounts,
    });
  };

  useEffect(() => {
    setState({
      lowFreqFilterType: configStore.lowFreqFilterType,
      minLocalCountsToShow: configStore.minLocalCountsToShow,
      minGlobalCountsToShow: configStore.minGlobalCountsToShow,
      maxGroupCounts: configStore.maxGroupCounts,
      hasChanged: false,
    });
  }, [
    configStore.lowFreqFilterType,
    configStore.minLocalCountsToShow,
    configStore.minGlobalCountsToShow,
    configStore.maxGroupCounts,
  ]);

  return (
    <SelectContainer>
      <SelectItem>
        <Radio>
          <input
            type="radio"
            value={LOW_FREQ_FILTER_TYPES.LOCAL_COUNTS}
            onChange={setFilterType}
            checked={
              state.lowFreqFilterType === LOW_FREQ_FILTER_TYPES.LOCAL_COUNTS
            }
          />
          <RadioLabel
            itemSelected={
              state.lowFreqFilterType === LOW_FREQ_FILTER_TYPES.LOCAL_COUNTS
            }
          >
            Minimum local counts
          </RadioLabel>
        </Radio>
        <Spacer />
        <NumberInput>
          <input
            type="number"
            value={state.minLocalCountsToShow}
            onChange={setMinLocalCountsToShow}
            disabled={
              state.lowFreqFilterType !== LOW_FREQ_FILTER_TYPES.LOCAL_COUNTS
            }
            min={1}
            step={1}
          />
        </NumberInput>
        <QuestionButton
          data-tip={`<p><b>${configStore.getGroupLabel()}s</b> with less than <b>${
            state.minLocalCountsToShow
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
              state.lowFreqFilterType === LOW_FREQ_FILTER_TYPES.GLOBAL_COUNTS
            }
          />
          <RadioLabel
            itemSelected={
              state.lowFreqFilterType === LOW_FREQ_FILTER_TYPES.GLOBAL_COUNTS
            }
          >
            Minimum global counts
          </RadioLabel>
        </Radio>
        <Spacer />
        <NumberInput>
          <input
            type="number"
            value={state.minGlobalCountsToShow}
            onChange={setMinGlobalCountsToShow}
            disabled={
              state.lowFreqFilterType !== LOW_FREQ_FILTER_TYPES.GLOBAL_COUNTS
            }
            min={1}
            step={1}
          />
        </NumberInput>
        <QuestionButton
          data-tip={`<p><b>${configStore.getGroupLabel()}s</b> with less than <b>${
            state.minGlobalCountsToShow
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
              state.lowFreqFilterType === LOW_FREQ_FILTER_TYPES.GROUP_COUNTS
            }
          />
          <RadioLabel
            itemSelected={
              state.lowFreqFilterType === LOW_FREQ_FILTER_TYPES.GROUP_COUNTS
            }
          >
            Show Top N {configStore.getGroupLabel()}s
          </RadioLabel>
        </Radio>
        <Spacer />
        <NumberInput>
          <input
            type="number"
            value={state.maxGroupCounts}
            onChange={setMaxGroupCounts}
            disabled={
              state.lowFreqFilterType !== LOW_FREQ_FILTER_TYPES.GROUP_COUNTS
            }
            min={1}
            step={1}
          />
        </NumberInput>
        <QuestionButton
          data-tip={`<p>Only show the top <b>${
            state.maxGroupCounts
          }</b> <b>${configStore.getGroupLabel()}s</b> by counts in the selected locations</p>`}
          data-html="true"
          data-for="tooltip-filter-sidebar"
        />
      </SelectItem>
      <ConfirmButton show={state.hasChanged} onClick={updateStore}>
        Confirm
      </ConfirmButton>
    </SelectContainer>
  );
});

export default FilterDataIntoOther;

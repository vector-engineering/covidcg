import React from 'react';
// import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import styled from 'styled-components';

const SelectContainer = styled.div`
  margin: 7px 13px 5px 13px;
`;

const GroupKeySelectForm = styled.form`
  label {
    display: flex;
    flex-direction: row;
    align-items: center;
    margin-bottom: 5px;

    span {
      flex-shrink: 0;
      margin-right: 10px;
    }

    select {
      padding: 1px 5px;
      flex-grow: 1;
      width: 100%;
      border-radius: 3px;
      margin: 0px;
    }
  }
`;

const RadioForm = styled.form`
  display: flex;
  flex-direction: row;
  align-items: center;

  span {
    flex-shrink: 0;
    margin-right: 10px;
  }

  .radio-row {
    display: flex;
    flex-direction: row;
    align-items: center;
    justify-content: flex-start;
  }

  .radio-item {
    margin-left: 0.65em;
    &:first-child {
      margin-left: 0px;
    }

    display: flex;
    flex-direction: row;
    align-items: center;

    input {
      margin: 0px;
      margin-right: 0.5em;
    }

    span.disabled-text {
      font-weight: normal;
      font-size: 0.8em;
      color: #888;
      flex-shrink: 1;
      white-space: normal;
      line-height: normal;
      margin-left: 5px;
    }
  }
`;

const GroupBySelect = observer(() => {
  const { covidStore } = useStores();

  let handleGroupKeyChange = (event) => {
    covidStore.changeGrouping(event.target.value, covidStore.dnaOrAa);
  };

  let handleDnaOrAaChange = (event) => {
    covidStore.changeGrouping(covidStore.groupKey, event.target.value);
  };

  // group by options
  let groupByOptions = [
    { label: 'Lineage', value: 'lineage' },
    { label: 'Clade', value: 'clade' },
    { label: 'SNV', value: 'snp' },
    // { label: 'SNP Signature', value: 'snp_sig' },
  ];

  // Create option elements
  let optionElements = [];
  groupByOptions.forEach((opt) => {
    optionElements.push(
      <option key={opt.value} value={opt.value}>
        {opt.label}
      </option>
    );
  });

  let aaDisabledMessage = '';
  let aaDisabled = false;
  if (
    covidStore.coordinateMode !== 'gene' &&
    covidStore.coordinateMode !== 'protein'
  ) {
    aaDisabledMessage = ' (only for gene/protein)';
    aaDisabled = true;
  } else if (covidStore.groupKey !== 'snp') {
    if (
      covidStore.coordinateMode === 'gene' &&
      covidStore.selectedGene.gene === 'All Genes'
    ) {
      aaDisabled = true;
      aaDisabledMessage = ' (please select one gene)';
    } else if (
      covidStore.coordinateMode === 'protein' &&
      covidStore.selectedProtein.protein === 'All Proteins'
    ) {
      aaDisabled = true;
      aaDisabledMessage = ' (please select one protein)';
    }
  }

  return (
    <SelectContainer>
      <GroupKeySelectForm>
        <label>
          <span>Group sequences by</span>
          <select value={covidStore.groupKey} onChange={handleGroupKeyChange}>
            {optionElements}
          </select>
        </label>
      </GroupKeySelectForm>
      <RadioForm>
        <span>Mutation format</span>
        <div className="radio-row">
          <div className="radio-item">
            <input
              type="radio"
              id="dnaChoice"
              name="dnaOrAa"
              value="dna"
              checked={covidStore.dnaOrAa === 'dna'}
              onChange={handleDnaOrAaChange}
            ></input>
            <label htmlFor="dnaChoice">NT</label>
          </div>
          <div className="radio-item">
            <input
              type="radio"
              id="aaChoice"
              name="dnaOrAa"
              value="aa"
              checked={covidStore.dnaOrAa === 'aa'}
              disabled={aaDisabled}
              onChange={handleDnaOrAaChange}
            ></input>
            <label htmlFor="aaChoice">AA</label>
            <span className="disabled-text">{aaDisabledMessage}</span>
          </div>
        </div>
      </RadioForm>
    </SelectContainer>
  );
});

GroupBySelect.propTypes = {};
GroupBySelect.defaultProps = {};

export default GroupBySelect;

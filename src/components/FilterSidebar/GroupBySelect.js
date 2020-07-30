import React from 'react';
// import PropTypes from 'prop-types';
import { observer } from 'mobx-react';
import { useStores } from '../../stores/connect';
import styled from 'styled-components';

const SelectContainer = styled.div`
  padding: 7px 13px 0px 13px;
  box-shadow: 0px 2px 4px #aaa;
  background-color: #fcfcfc;
`;

const RadioForm = styled.form`
  display: flex;
  flex-direction: column;
  align-items: stretch;

  margin-bottom: 10px;
  // margin-bottom: 5px;

  // background-color: #fff;
  // padding: 5px 10px;
  // border: 1px solid #ccc;

  span {
    flex-shrink: 0;
    margin-right: 10px;
  }

  .radio-row {
    display: flex;
    flex-direction: row;
    align-items: center;
    justify-content: flex-start;

    margin-top: 3px;
    margin-left: 12px;
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
      <RadioForm>
        <span className="form-title">Group sequences by</span>
        <div className="radio-row">
          <div className="radio-item">
            <label>
              <input
                className="radio-input"
                type="radio"
                value="lineage"
                checked={covidStore.groupKey === 'lineage'}
                onChange={handleGroupKeyChange}
              />
              <span>Lineage</span>
            </label>
          </div>
          <div className="radio-item">
            <label>
              <input
                className="radio-input"
                type="radio"
                value="clade"
                checked={covidStore.groupKey === 'clade'}
                onChange={handleGroupKeyChange}
              />
              <span>Clade</span>
            </label>
          </div>
          <div className="radio-item">
            <label>
              <input
                className="radio-input"
                type="radio"
                value="snp"
                checked={covidStore.groupKey === 'snp'}
                onChange={handleGroupKeyChange}
              />
              <span>SNV</span>
            </label>
          </div>
        </div>
      </RadioForm>
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

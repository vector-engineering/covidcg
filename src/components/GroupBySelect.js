import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';

const SelectContainer = styled.div`
  margin: 5px 5px 0px 5px;
  padding: 5px 8px;
`;

const GroupKeySelectForm = styled.form`
  display: flex;
  flex-direction: column;
  padding-right: 20px;
  label {
    display: flex;
    flex-direction: column;
    align-items: flex-start;
    justify-content: flex-start;
  }
  select {
    padding: 1px 5px;
    margin-left: 0.75em;
    flex-grow: 1;
    margin-top: 5px;
    width: 100%;
    border-radius: 3px;
  }
`;

const RadioForm = styled.form`
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  justify-content: flex-start;
  margin-top: 0.85em;

  p {
    // margin-right: 0.1em;
    margin: 0px;
  }

  .radio-row {
    flex-grow: 1;
    display: flex;
    flex-direction: row;
    align-items: center;
    justify-content: flex-start;
    margin-top: 5px;
  }

  .radio-item {
    margin-left: 0.65em;
    display: flex;
    flex-direction: row;
    align-items: center;

    input {
      margin: 0px;
      margin-right: 0.5em;
    }
  }
`;

const GroupBySelect = ({ groupKey, dnaOrAa, onChange }) => {
  let handleGroupKeyChange = (event) => {
    onChange(event.target.value, dnaOrAa);
  };

  let handleDnaOrAaChange = (event) => {
    onChange(groupKey, event.target.value);
  };

  // group by options
  let groupByOptions = [
    { label: 'Lineage', value: 'lineage' },
    { label: 'SNP', value: 'snp' },
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

  return (
    <SelectContainer>
      <GroupKeySelectForm>
        <label>
          Group sequences by:
          <select value={groupKey} onChange={handleGroupKeyChange}>
            {optionElements}
          </select>
        </label>
      </GroupKeySelectForm>
      <RadioForm>
        <p>Mutation format:</p>
        <div className="radio-row">
          <div className="radio-item">
            <input
              type="radio"
              id="dnaChoice"
              name="dnaOrAa"
              value="dna"
              checked={dnaOrAa === 'dna'}
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
              checked={dnaOrAa === 'aa'}
              onChange={handleDnaOrAaChange}
            ></input>
            <label htmlFor="aaChoice">AA</label>
          </div>
        </div>
      </RadioForm>
    </SelectContainer>
  );
};

GroupBySelect.propTypes = {
  groupKey: PropTypes.string,
  dnaOrAa: PropTypes.string,
  onChange: PropTypes.func,
};

GroupBySelect.defaultProps = {
  groupKey: 'lineage',
  dnaOrAa: 'dna',
  onChange: () => {},
};

export default GroupBySelect;

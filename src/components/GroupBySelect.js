import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';

const SelectContainer = styled.div`
  margin: 5px;
  padding: 5px 8px;
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
      <form>
        <label>
          Select a grouping key
          <select value={groupKey} onChange={handleGroupKeyChange}>
            {optionElements}
          </select>
        </label>
      </form>
      <form>
        <p>Select Mutation Format</p>
        <div>
          <input
            type="radio"
            id="dnaChoice"
            name="dnaOrAa"
            value="dna"
            checked={dnaOrAa === 'dna'}
            onChange={handleDnaOrAaChange}
          ></input>
          <label htmlFor="dnaChoice">NT</label>
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
      </form>
    </SelectContainer>
  );
};

GroupBySelect.propTypes = {
  groupKey: PropTypes.string.isRequired,
  dnaOrAa: PropTypes.string.isRequired,
  onChange: PropTypes.func,
};

GroupBySelect.defaultProps = {
  groupKey: 'lineage',
  dnaOrAa: 'dna',
  onChange: () => {},
};

export default GroupBySelect;

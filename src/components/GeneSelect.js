import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';

const SelectContainer = styled.div`
  margin: 5px;
  padding: 5px 8px;
  margin-top: 0px;
  padding-top: 0px;
`;
const GeneSelectForm = styled.form`
  font-weight: 700;

  label {
    display: flex;
    flex-direction: row;
    align-items: center;
    justify-content: flex-start;
  }

  select {
    background-color: white;
    flex-grow: 1;
    margin-left: 0.65em;
    padding: 1px 5px;
  }
`;
const PositionContainer = styled.div`
  display: grid;
  grid-template-columns: [col1] 50% [col2] 50% [col3];
  grid-template-rows: [row1] auto [row2];
  font-weight: 500;
  font-size: 0.9em;
  margin-top: 3px;
`;
const PosFrom = styled.div`
  grid-column: col1/col2;
  grid-row: row1/row2;
  padding-right: 8px;
  border-right: 1px solid #aaa;
`;
const PosTo = styled.div`
  grid-column: col2/col3;
  grid-row: row1/row2;
  padding-left: 8px;
`;

const GeneSelect = ({ genes, selectedGene, onChange }) => {
  const handleChange = (event) => {
    onChange(event.target.value);
  };

  // Create option elements
  let optionElements = [];
  genes.forEach((opt) => {
    optionElements.push(
      <option key={opt.value} value={opt.value}>
        {opt.label}
      </option>
    );
  });

  return (
    <SelectContainer>
      <GeneSelectForm>
        <label>
          Gene:
          <select value={selectedGene.gene} onChange={handleChange}>
            {optionElements}
          </select>
        </label>
      </GeneSelectForm>
      <PositionContainer>
        <PosFrom>From: {selectedGene.start}</PosFrom>
        <PosTo>To: {selectedGene.end}</PosTo>
      </PositionContainer>
    </SelectContainer>
  );
};

GeneSelect.propTypes = {
  genes: PropTypes.array.isRequired,
  selectedGene: PropTypes.object.isRequired,
  onChange: PropTypes.func,
};

GeneSelect.defaultProps = {
  onChange: () => {},
};

export default GeneSelect;

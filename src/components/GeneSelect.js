import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';

const SelectContainer = styled.div`
  margin: 5px;
  padding: 5px 8px;
`;
const GeneSelectLabel = styled.label`
  font-weight: 700;
`;
const GeneSelectSelect = styled.select`
  width: 100%;
  background-color: white;
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
      <form>
        <GeneSelectLabel>
          Select a gene to analyze
          <GeneSelectSelect value={selectedGene.gene} onChange={handleChange}>
            {optionElements}
          </GeneSelectSelect>
        </GeneSelectLabel>
      </form>
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

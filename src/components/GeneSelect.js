import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';

const SelectContainer = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;
  justify-content: flex-start;

  margin: 5px;
  padding: 5px 8px;
  margin-top: 0px;
  padding-top: 0px;
  margin-bottom: 5px;
`;
const GeneSelectForm = styled.form`
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  padding-right: 20px;

  label {
    display: flex;
    flex-direction: row;
    align-items: center;
    justify-content: flex-start;
  }

  select {
    margin-top: 5px;
    background-color: white;
    flex-grow: 1;
    margin-left: 0.65em;
    padding: 1px 5px;
    width: 100%;
  }
`;
const PositionContainer = styled.div`
  display: flex;
  font-weight: normal;
  font-size: 0.9em;
  margin-top: 3px;
  margin-left: 8px;
  margin-right: 20px;
  justify-content: space-between;
`;
const PosFrom = styled.div`
  padding-right: 8px;
  border-right: 1px solid #aaa;
  flex-grow: 1;
`;
const PosTo = styled.div`
  padding-left: 8px;
  flex-grow: 1;
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
        <label>Gene:</label>
        <select value={selectedGene.gene} onChange={handleChange}>
          {optionElements}
        </select>
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

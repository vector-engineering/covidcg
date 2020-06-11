import React from 'react';
import PropTypes from 'prop-types';
import styled from 'styled-components';

const LetterCellDiv = styled.div`
  font-family: monospace;
  font-weight: 500;
  font-size: 1.25em;
  width: 25px;
  height: 100%;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
`;

const LetterCell = ({ value, bgColor }) => {
  return (
    <LetterCellDiv style={{ backgroundColor: bgColor }}>{value}</LetterCellDiv>
  );
};

LetterCell.propTypes = {
  value: PropTypes.string.isRequired,
  bgColor: PropTypes.string,
};
LetterCell.defaultProps = {
  bgColor: 'transparent',
};

export default LetterCell;

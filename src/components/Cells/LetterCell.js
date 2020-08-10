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

  background-color: ${({ bgColor }) => bgColor};
  color: ${({ textColor }) => textColor};
`;

const LetterCell = ({ value, bgColor, textColor }) => {
  return (
    <LetterCellDiv bgColor={bgColor} textColor={textColor}>
      {value}
    </LetterCellDiv>
  );
};

LetterCell.propTypes = {
  value: PropTypes.string,
  bgColor: PropTypes.string,
  textColor: PropTypes.string,
};
LetterCell.defaultProps = {
  value: '',
  bgColor: 'transparent',
  textColor: '#000',
};

export default LetterCell;

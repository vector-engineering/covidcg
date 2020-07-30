import React from 'react';
import styled from 'styled-components';

const QuestionButtonContainer = ({ ...props }) => {
  return <span {...props}>?</span>;
};

const QuestionButton = styled(QuestionButtonContainer)`
  font-family: monospace;
  font-size: 1em;
  line-height: normal;

  margin-left: 8px;
  padding: 2px 5px;
  background-color: #fff;
  border: 1px solid #ccc;
  border-radius: 2px;
  &:hover {
    background-color: #f8f8f8;
  }
`;

export default QuestionButton;

import React from 'react';
import styled from 'styled-components';

const EmptyDataTableContainer = styled.div`
  min-height: 300px;
  width: 100%;
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: center;

  span {
    font-weight: normal;
    font-size: 1.5em;
  }
`;

const EmptyDataTable = () => {
  return (
    <EmptyDataTableContainer>
      <span>No sequences selected</span>
    </EmptyDataTableContainer>
  );
};

export default EmptyDataTable;

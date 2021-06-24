import styled from 'styled-components';

export const VOCGridTitle = styled.span`
  font-size: 1em;
  grid-row: 1;
  justify-self: center;
  font-weight: bold;
`;

export const GridItem = styled.div`
  padding-left: 10px;
`;

export const VOCListContainer = styled.div`
  display: flex;
  flex-direction: column;
  padding: 10px;
  overflow: hidden;
  font-weight: normal;
`;

export const VOCItemGrid = styled.div`
  display: grid;
  grid-template-columns: repeat(2, 1fr);
  grid-template-rows: auto;
  grid-column-gap: 10px;
`;

export const VOCItemContainer = styled.div`
  position: relative;

  input {
    vertical-align: 'baseline';
  }
`;

export const VOCItemName = styled.span`
  &:hover {
    cursor: pointer;
  }
`;

export const DropdownGrid = styled.div`
  display: grid;
  top: 100%;
  z-index: 2;
  grid-template-columns: 180px;
  grid-template-rows: auto;
  background-color: #ffffff;
`;

export const DropdownLabel = styled.span`
  grid-column: 1;
  grid-row: auto;
`;

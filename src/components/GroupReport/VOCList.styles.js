import styled from 'styled-components';

const StyledLi = styled.li`
  float: left;
  height: 100%;
  padding: 10px;
`;

const VOCGridTitle = styled.span`
  font-size: 16px;
  grid-row: 1;
  justify-self: center;
`;

const GridItem = styled(StyledLi)`
  grid-row: auto;
  list-style-type: none;
  justify-self: left;
  position: relative;
`;

export const VOCGridItem = styled(GridItem)`
  grid-column: 1;
`;

export const VOIGridItem = styled(GridItem)`
  grid-column: 2;
`;

export const VOCListContainer = styled.div`
  display: flex;
  flex-direction: column;
  padding: 10px;
  overflow: hidden;
`;

export const VOCListHeader = styled.div`
  text-align: center;
`;

export const VOCListTitle = styled.div`
  font-size: 18px;
`;

export const VOCItemGrid = styled.div`
  display: grid;
  grid-template-columns: repeat(2, 180px);
  grid-template-rows: auto;
`;

export const VOCTitle = styled(VOCGridTitle)`
  grid-column: 1;
`;

export const VOITitle = styled(VOCGridTitle)`
  grid-column: 2;
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

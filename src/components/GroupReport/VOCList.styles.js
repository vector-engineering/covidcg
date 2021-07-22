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
  display: grid;
  grid-template-columns: 20px 100px;
  grid-gap: 5px
  align-items: center;
  background-color: ${({ selected }) => (selected ? '#117733' : '#DDDDDD')};
  margin-top: 5px;
  padding: 5px;

  position: relative;

  &:hover {
    cursor: pointer;
  }
`;

export const VOCItemName = styled.span`
  float: left;
  color: ${({ selected }) => (selected ? 'white' : 'black')};
`;

export const VOCBadgeContainer = styled.div`
  display: grid;
  margin: auto;
  max-width: 11px;
  max-height: 11px;
  border: 1px solid black;
  background-color: black;
  grid-template-columns: repeat(2, 5px);
  grid-template-rows: repeat(2, 5px);
  grid-gap: 1px;
`;

export const VOCBadge = styled.div`
  background-color: ${({ color }) => color};
  grid-row: ${({ row }) => row};
  grid-column: ${({ col }) => col};
`;
VOCBadge.defaultProps = { color: 'white', row: 1, col: 1 };

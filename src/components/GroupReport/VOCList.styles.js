import styled from 'styled-components';

export const VOCGridTitle = styled.span`
  font-size: 1em;
  grid-row: 1;
  font-weight: bold;
  text-align: center;
  grid-column: span ${({ colSpan }) => colSpan};
`;
VOCGridTitle.defaultProps = {
  colSpan: 1,
};

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
  grid-template-columns: repeat(4, 1fr);
  grid-template-rows: auto;
  grid-column-gap: 10px;
  justify-items: center;
`;

export const VOCItemContainer = styled.div`
  display: grid;
  grid-template-columns: 11px 1fr;
  grid-gap: 5px
  align-items: center;
  background-color: ${({ selected }) => (selected ? '#009988' : '#eeeeee')};
  margin-top: 5px;
  padding: 1px 5px;
  font-size: 0.85em;
  border-radius: 5px;

  position: relative;

  &:hover {
    cursor: pointer;
  }
`;

export const VOCItemName = styled.span`
  margin-left: 5px;
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

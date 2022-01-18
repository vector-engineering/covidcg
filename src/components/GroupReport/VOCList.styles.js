import styled, { css } from 'styled-components';

export const VOCGridTitle = styled.span`
  font-size: 1.5em;
`;

export const VOCTableContainer = styled.div`
  display: flex;
  flex-direction: column;
  width: 100%;
`;

export const VOCTableHeader = styled.div`
  display: flex;
  flex-direciton: column;
`;

export const VOCTableToggle = styled.div`
  transform: rotate(0deg);
  margin: auto auto auto 10px;
  width: 0;
  height: 0;
  border-left: 5px solid transparent;
  border-right: 5px solid transparent;
  border-top: 10px solid #444444;
  cursor: pointer;

  ${(props) =>
    props.expanded &&
    css`
      transition-timing-function: ease-in-out;
      transform: rotate(180deg);
    `};

  &:hover {
    border-top: 10px solid #aaaaaa;
  }
`;

export const VOCTableRow = styled.div`
  display: flex;
  flex-direction: column;
  margin-top: 0.5em;
`;

export const VOCTableContent = styled.div`
  display: flex;
  flex-direction: row;
  flex-wrap: wrap;
`;

export const VOCItemContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  background-color: ${({ selected }) => (selected ? '#009988' : '#eeeeee')};
  margin: 5px;
  padding: 5px 5px;
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

  &:after {
    margin-left: 5px;
    content: ${({ whoLabel }) => (whoLabel ? JSON.stringify(whoLabel) : '')};
  }
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

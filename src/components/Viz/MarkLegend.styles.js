import styled from 'styled-components';

export const MarkLegendContainer = styled.div`
  display: flex;
  flex-direction: column;
  align-items: stretch;
`;

export const LegendTitle = styled.div``;

export const MarkItemList = styled.div`
  display: grid;
  grid-template-columns: repeat(3, auto);
  grid-column-gap: 3px;
  grid-auto-flow: row;
`;

export const MarkItemContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
`;

export const MarkColorBox = styled.div`
  width: 10px;
  height: 10px;
  background-color: ${({ color }) => color};
  margin-right: 4px;
`;
MarkColorBox.defaultProps = {
  color: '#AAA',
};

export const MarkLabel = styled.div`
  font-size: 0.85em;
  line-height: normal;
  white-space: nowrap;
`;
